"""Explore fixed ortho-site contact descriptors; total forces are not pair forces."""
from __future__ import annotations
import itertools,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from scipy.stats import pearsonr,spearmanr

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.extend([str(ROOT/'base_analysis/src'),str(ROOT/'analysis2_force/src')])
from force_motion import _trajectory,_align_to_ts
from reporting import markdown
from plotting4 import SYSTEM_COLORS
import matplotlib.pyplot as plt

KCAL=627.509474;BOHR_A=.529177210903

def angle(x,y,z):
    a=x-y;b=z-y
    return np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',a,b)/(np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1)),-1,1)))


def ortho_substituent(xyz,z,carbon):
    candidates=[i for i in range(len(z)) if z[i] in (1,9) and np.linalg.norm(xyz[i]-xyz[carbon-1])<1.6]
    if len(candidates)!=1:raise ValueError(f'Nonunique substituent at C{carbon}: {candidates}')
    return candidates[0]


def pair_projection(xyz,forces,a,b,velocity):
    delta=xyz[:,b]-xyz[:,a];q=np.linalg.norm(delta,axis=1);u=delta/q[:,None]
    projection=.5*np.einsum('ij,ij->i',forces[:,b]-forces[:,a],u)
    dq_bohr=np.einsum('ij,ij->i',velocity[:,b]-velocity[:,a],u)
    return q,projection,projection*dq_bohr


def calculate_profiles(system,frame,cfg):
    keys,atoms,xyz,forces=_trajectory(frame);r=keys.signed_irc_coordinate.to_numpy(float)
    ts=int(np.argmin(abs(r)));z=frame[frame.signed_irc_coordinate==r[ts]].sort_values('atom_number').atomic_number.to_numpy(int)
    n=44 if system=='MoriPF' else 42
    hs=[i for i in range(len(z)) if z[i]==1 and np.linalg.norm(xyz[ts,i]-xyz[ts,n])<1.3]
    assert len(hs)==2
    near=ortho_substituent(xyz[ts],z,cfg['near_ring_C']);opp=ortho_substituent(xyz[ts],z,cfg['opposite_ring_C']);o=cfg['opposite_O']-1
    ordered=sorted(hs,key=lambda h:np.linalg.norm(xyz[ts,near]-xyz[ts,h]))
    xyz,forces=_align_to_ts(xyz,forces,ts)
    velocity=np.gradient(xyz,r,axis=0,edge_order=2)/BOHR_A
    atomic_C=np.einsum('pij,pij->pi',forces,velocity);norm=np.linalg.norm(forces,axis=2)
    symbols={1:'H',9:'F'}
    audit=dict(system=system,near_ring_C=cfg['near_ring_C'],near_atom=near+1,near_element=symbols[z[near]],primary_amine_H=ordered[0]+1,
        alternative_amine_H=ordered[1]+1,opposite_ring_C=cfg['opposite_ring_C'],opposite_atom=opp+1,opposite_element=symbols[z[opp]],opposite_O=o+1)
    qopp,Qopp,Aopp=pair_projection(xyz,forces,opp,o,velocity)
    frames=[]
    for choice,h in zip(['primary','alternative_amine_H'],ordered):
        qnear,Qnear,Anear=pair_projection(xyz,forces,near,h,velocity)
        subset=sorted({near,h,opp,o});assert len(subset)==4
        energies=frame.groupby('signed_irc_coordinate').energy_hartree.first().reindex(r).to_numpy()
        df=pd.DataFrame(dict(system=system,H_choice=choice,R=r,E_hartree=energies,CN_A=np.linalg.norm(xyz[:,19]-xyz[:,n],axis=1),
            near_atom=near+1,amine_H=h+1,opposite_atom=opp+1,opposite_O=o+1,
            near_distance_A=qnear,opposite_distance_A=qopp,
            near_NHX_angle_deg=angle(xyz[:,n],xyz[:,h],xyz[:,near]),
            opposite_CXO_angle_deg=angle(xyz[:,cfg['opposite_ring_C']-1],xyz[:,opp],xyz[:,o]),
            near_projected_force_hartree_bohr=Qnear,opposite_projected_force_hartree_bohr=Qopp,
            near_force_motion_hartree_R=Anear,opposite_force_motion_hartree_R=Aopp,
            near_X_force_magnitude=norm[:,near],amine_H_force_magnitude=norm[:,h],opposite_X_force_magnitude=norm[:,opp],opposite_O_force_magnitude=norm[:,o],
            near_atomic_force_magnitude_sum=norm[:,[near,h]].sum(axis=1),
            opposite_atomic_force_magnitude_sum=norm[:,[opp,o]].sum(axis=1),
            selected_atomic_force_magnitude_sum=norm[:,subset].sum(axis=1),all_atomic_force_magnitude_sum=norm.sum(axis=1),
            near_X_atomic_C_hartree_R=atomic_C[:,near],amine_H_atomic_C_hartree_R=atomic_C[:,h],
            opposite_X_atomic_C_hartree_R=atomic_C[:,opp],opposite_O_atomic_C_hartree_R=atomic_C[:,o],
            near_atomic_C_hartree_R=atomic_C[:,[near,h]].sum(axis=1),
            opposite_atomic_C_hartree_R=atomic_C[:,[opp,o]].sum(axis=1),
            remaining_atomic_C_hartree_R=atomic_C[:,[i for i in range(len(z)) if i not in subset]].sum(axis=1),
            selected_atomic_C_hartree_R=atomic_C[:,subset].sum(axis=1),all_atomic_C_hartree_R=atomic_C.sum(axis=1)))
        frames.append(df)
    return pd.concat(frames,ignore_index=True),audit


def window_metrics(profile,cfg,energy):
    results=[]
    for (system,choice),g in profile.groupby(['system','H_choice'],sort=False):
        g=g.sort_values('R');ts=g.loc[g.R.abs().idxmin()]
        pre=g[g.R<0].dropna(subset=['near_projected_force_hartree_bohr','opposite_projected_force_hartree_bohr'])
        assert (np.diff(pre.CN_A)<0).all(),f'{system}: pre-TS CN is not monotonic'
        for window,(start,end) in cfg['windows'].items():
            if not (pre.CN_A.max()>=start>end>=pre.CN_A.min()):raise ValueError(f'{system}: common window not covered')
            # R grid includes interpolated boundaries plus every computed interior point.
            r0=np.interp(start,pre.CN_A.to_numpy()[::-1],pre.R.to_numpy()[::-1]);r1=np.interp(end,pre.CN_A.to_numpy()[::-1],pre.R.to_numpy()[::-1])
            rr=np.r_[r0,pre.R[(pre.R>r0)&(pre.R<r1)],r1]
            def onr(col):return np.interp(rr,pre.R,pre[col])
            def integral(col):return float(np.trapezoid(onr(col),rr))
            # Equal progress weighting for force averages and shares, independent of IRC point spacing.
            cn=np.linspace(end,start,501)
            def oncn(col):return np.interp(cn,pre.CN_A.to_numpy()[::-1],pre[col].to_numpy()[::-1])
            def average(col,absolute=False):
                v=oncn(col);return float(np.trapezoid(abs(v) if absolute else v,cn)/(start-end))
            nearabs=average('near_projected_force_hartree_bohr',True);oppabs=average('opposite_projected_force_hartree_bohr',True)
            nearatomic=average('near_atomic_force_magnitude_sum');oppatomic=average('opposite_atomic_force_magnitude_sum')
            selectedatom=average('selected_atomic_force_magnitude_sum');allatom=average('all_atomic_force_magnitude_sum')
            costnear=-integral('near_force_motion_hartree_R')*KCAL;costopp=-integral('opposite_force_motion_hartree_R')*KCAL
            e0=float(np.interp(r0,g.R,g.E_hartree));e1=float(np.interp(r1,g.R,g.E_hartree))
            results.append(dict(
                entry_energy_kcal_mol=(e0-energy.loc[system,'separated_reactants_E_hartree'])*KCAL,
                window_energy_rise_kcal_mol=(e1-e0)*KCAL,
                after_window_to_TS_kcal_mol=(energy.loc[system,'TS_E_hartree']-e1)*KCAL,
                work_energy_residual_kcal_mol=-integral('all_atomic_C_hartree_R')*KCAL-(e1-e0)*KCAL,
                system=system,window=window,H_choice=choice,R_start=r0,R_end=r1,CN_start_A=start,CN_end_A=end,
                near_X_mean_force_mEh_bohr=average('near_X_force_magnitude')*1000,amine_H_mean_force_mEh_bohr=average('amine_H_force_magnitude')*1000,
                opposite_X_mean_force_mEh_bohr=average('opposite_X_force_magnitude')*1000,opposite_O_mean_force_mEh_bohr=average('opposite_O_force_magnitude')*1000,
                near_mean_abs_projected_force_mEh_bohr=nearabs*1000,opposite_mean_abs_projected_force_mEh_bohr=oppabs*1000,
                near_mean_signed_projected_force_mEh_bohr=average('near_projected_force_hartree_bohr')*1000,
                opposite_mean_signed_projected_force_mEh_bohr=average('opposite_projected_force_hartree_bohr')*1000,
                opposite_projected_magnitude_share=oppabs/(nearabs+oppabs),
                opposite_atomic_force_share=oppatomic/(nearatomic+oppatomic),selected_atomic_force_share=selectedatom/allatom,
                near_pair_diagnostic_work_kcal_mol=costnear,opposite_pair_diagnostic_work_kcal_mol=costopp,
                pair_diagnostic_work_kcal_mol=costnear+costopp,
                **{name+'_atomic_work_kcal_mol':-integral(name+'_atomic_C_hartree_R')*KCAL for name in ['near_X','amine_H','opposite_X','opposite_O','near','opposite','remaining']},
                selected_atomic_work_kcal_mol=-integral('selected_atomic_C_hartree_R')*KCAL,
                all_atomic_work_kcal_mol=-integral('all_atomic_C_hartree_R')*KCAL,
                near_TS_distance_A=ts.near_distance_A,opposite_TS_distance_A=ts.opposite_distance_A,opposite_TS_CXO_angle_deg=ts.opposite_CXO_angle_deg,
                absolute_electronic_barrier_kcal_mol=energy.loc[system,'absolute_electronic_barrier_kcal_mol'],
                absolute_RRHO_G_barrier_kcal_mol=energy.loc[system,'absolute_RRHO_G_barrier_kcal_mol']))
    metrics=pd.DataFrame(results)
    return metrics


def correlations(metrics,cfg):
    rows=[]
    for (window,choice),g in metrics.groupby(['window','H_choice'],sort=False):
        for metric,target in itertools.product(cfg['correlation_descriptors'],['absolute_electronic_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol']):
            x=g[metric].to_numpy();y=g[target].to_numpy()
            r=float(pearsonr(x,y).statistic);rho=float(spearmanr(x,y).statistic)
            # Fit the uncentered descriptor and barrier with an intercept.
            loo=[];rleave=[]
            for i in range(4):
                mask=np.arange(4)!=i
                slope,intercept=np.polyfit(x[mask],y[mask],1)
                loo.append(slope*x[i]+intercept)
                rleave.append(float(pearsonr(x[mask],y[mask]).statistic))
            rows.append(dict(window=window,H_choice=choice,descriptor=metric,target=target,n_systems=4,
                pearson_r=r,spearman_r=rho,LOO_prediction_RMSE_kcal_mol=float(np.sqrt(np.mean((np.array(loo)-y)**2))),
                leave_one_out_r_min=min(rleave),leave_one_out_r_max=max(rleave),
                exploratory_only=True,independent_validation=False))
    return pd.DataFrame(rows)


def main():
    cfg=yaml.safe_load((HERE/'config/two_sided_contacts.yaml').read_text());out=HERE/'outputs';csv=out/'csv';fig=out/'figures';reports=out/'reports'
    energy=pd.read_csv(csv/'energy_comparison.csv').set_index('system');pieces=[];audits=[]
    for system in energy.index:
        profile,audit=calculate_profiles(system,pd.read_csv(csv/f'parsed_oriented_{system}.csv.gz'),cfg)
        pieces.append(profile);audits.append(audit)
    profiles=pd.concat(pieces,ignore_index=True);metrics=window_metrics(profiles,cfg,energy);corr=correlations(metrics,cfg)
    pd.DataFrame(audits).to_csv(csv/'two_sided_atom_mapping.csv',index=False)
    profiles.to_csv(csv/'two_sided_force_profiles.csv',index=False);metrics.to_csv(csv/'two_sided_window_metrics.csv',index=False);corr.to_csv(csv/'two_sided_correlations.csv',index=False)
    primary=metrics[(metrics.window==cfg['primary_window'])&(metrics.H_choice=='primary')]
    print(primary[['system','near_mean_abs_projected_force_mEh_bohr','opposite_mean_abs_projected_force_mEh_bohr','opposite_projected_magnitude_share','opposite_atomic_force_share','selected_atomic_force_share','selected_atomic_work_kcal_mol','pair_diagnostic_work_kcal_mol','opposite_TS_distance_A','opposite_TS_CXO_angle_deg']].to_string(index=False))
    print(corr[(corr.window==cfg['primary_window'])&(corr.H_choice=='primary')].to_string(index=False))
    from two_sided_report import create_report
    create_report(profiles,metrics,corr,pd.DataFrame(audits),cfg,fig,reports)

if __name__=='__main__':main()
