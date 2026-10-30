"""Durable continuation with at most two systems in parallel.

No registration, additional software, optimization, or all-IRC expansion.
An existing RAFA14 worker occupies one slot while another system runs.
"""
import argparse
import ctypes
import csv
import json
import os
from pathlib import Path
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from irc_eda.scripts.smd_distortion_pilot import OUT, REPORT, save, runjob, validate, analyze


def wait_pid(pid):
    if not pid:return
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.OpenProcess.argtypes=[ctypes.c_ulong,ctypes.c_int,ctypes.c_ulong]
    kernel.OpenProcess.restype=ctypes.c_void_p
    kernel.WaitForSingleObject.argtypes=[ctypes.c_void_p,ctypes.c_ulong]
    kernel.CloseHandle.argtypes=[ctypes.c_void_p]
    handle=kernel.OpenProcess(0x00100000,False,pid)
    if not handle:
        if ctypes.get_last_error() not in (0,87):raise ctypes.WinError(ctypes.get_last_error())
        return
    try:
        while True:
            result=kernel.WaitForSingleObject(handle,1000)
            if result==0:return
            if result!=258:raise ctypes.WinError(ctypes.get_last_error())
    finally:kernel.CloseHandle(handle)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--wait-pid',type=int);parser.add_argument('--gaussian-pid',type=int);args=parser.parse_args()
    status=OUT/'workflow_status.json'
    lock=OUT/'workflow.lock'
    with lock.open('x') as f:f.write(str(os.getpid()))
    record={'analysis':'fixed-fragment SMD distortion/interaction differences','status':'WAITING_FOR_EXISTING_RAFA14_WORKER','started_unix':time.time(),'systems':['RAFA14','MoriPF','RAFA15','RAFA16'],'completed_systems':[],'scope':'3 points per system; 24 fragment SPs maximum, successful attempts reused','wait_pid':args.wait_pid,'gaussian_pid':args.gaussian_pid,'max_parallel_systems':2,'cores_per_job':4,'memory_gb_per_job':4,'active_jobs':{}}
    save(status,record)
    print('Waiting for existing RAFA14 worker',args.wait_pid,flush=True)
    mutex=threading.Lock()
    stopped=threading.Event()
    def run_system(system):
        try:
            if system=='RAFA14':
                with mutex:
                    record['active_jobs'][system]={'existing_worker_pid':args.wait_pid,'state':'RUNNING_EXTERNALLY'}
                    save(status,record)
                # Reserve this slot until the original worker and its children
                # finish. Never launch another job in its folders concurrently.
                if args.wait_pid:wait_pid(args.wait_pid)
                else:wait_pid(args.gaussian_pid)
            manifest=json.loads((OUT/system/'manifest.json').read_text())
            for job in manifest['jobs']:
                if stopped.is_set():raise RuntimeError('Stopped after another system failed')
                with mutex:
                    record.update(status='RUNNING')
                    record['active_jobs'][system]=job['folder'];save(status,record)
                folder=ROOT/job['folder'];validation=folder/'validation.json'
                inp=(folder/'fragment.gjf').read_text()
                assert '%nprocshared=4' in inp.lower()
                assert any('%mem='+size in inp.lower() for size in ('2gb','4gb'))
                # If a previous Python worker exited after Gaussian finished,
                # validate that output rather than overwrite it or recompute.
                if (folder/'fragment.out').exists():
                    existing=json.loads(validation.read_text()) if validation.exists() else {}
                    if existing.get('status')!='SUCCESS':
                        recovered=validate(job)
                        recovered['recovered_after_worker_exit']=True;save(validation,recovered)
                runjob(job)
            # Plotting and shared status updates are serialized.
            with mutex:
                analyze(system)
                record['completed_systems'].append(system)
                record['active_jobs'].pop(system,None);save(status,record)
            print('REPORTED',system,str(REPORT/system/'report_ja.md'),flush=True)
        except Exception:
            stopped.set()
            raise
    try:
        pending_systems=iter(record['systems'])
        with ThreadPoolExecutor(max_workers=2) as pool:
            pending={pool.submit(run_system,next(pending_systems)) for _ in range(2)}
            while pending:
                done,pending=wait(pending,return_when=FIRST_COMPLETED)
                for future in done:future.result()
                for _ in done:
                    system=next(pending_systems,None)
                    if system is not None and not stopped.is_set():
                        pending.add(pool.submit(run_system,system))
        rows=[]
        for system in record['systems']:
            for row in json.loads((REPORT/system/'pilot.json').read_text()):rows.append({'system':system,**row})
        with (REPORT/'all_pilots.csv').open('w',newline='',encoding='utf8') as f:
            writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
        record.update(status='COMPLETE',finished_unix=time.time(),current_job=None,current_system=None)
        save(status,record)
    except Exception as exc:
        record.update(status='FAILED_STOPPED',error=repr(exc),stopped_unix=time.time());save(status,record)
        raise
    finally:
        lock.unlink()


if __name__=='__main__':main()
