"""Reproducible four-system analysis; run from any working directory."""
from __future__ import annotations
import copy, json, sys, platform
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import yaml
from scipy.stats import pearsonr,spearmanr

ROOT=Path(__file__).resolve().parents[1]; HERE=Path(__file__).resolve().parent
sys.path.extend([str(ROOT/'base_analysis/src'),str(ROOT/'analysis2_force/src')])
from inputs import restore_inputs,parse_stationary,serializable,prepare_irc,aligned,bond,sha256
from io_parsers import parse_irc_dipoles
from numerics4 import reaction_force_profile
from geometry_context import analyze_geometry_context,analyze_atomic_forces
from force_motion import analyze_force_motion
import plotting4 as plots
from energetics import energy_terms, validate_reference

KCAL=627.509474

def table(rows,path):
    result=pd.DataFrame(rows); result.to_csv(path,index=False);return result


def angle(x,y,z):
    a=x-y;b=z-y
    return float(np.degrees(np.arccos(np.clip(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b),-1,1))))


def mapping_from_ts(system,ts,n):
    xyz=ts['xyz'];zs=ts['atomic_numbers'];f=list(np.flatnonzero(zs==9)+1)
    hs=[int(i+1) for i in np.flatnonzero(zs==1) if np.linalg.norm(xyz[i]-xyz[n-1])<1.3]
    assert len(hs)==2
    # Select once at the TS; retain identity along the entire IRC.
    distance,ff,hh=min((np.linalg.norm(xyz[a-1]-xyz[b-1]),int(a),b) for a in f for b in hs)
    oh=min(hs,key=lambda h:np.linalg.norm(xyz[h-1]-xyz[18]))
    interference=int(min(f,key=lambda a:np.linalg.norm(xyz[a-1]-xyz[20])))
    mapping=dict(reaction_C=20,reaction_N=n,hb_substituent_F=ff,hb_substituent_H=hh,hb_coumarin_O=19,hb_coumarin_H=oh,interference_F=interference,interference_O=21)
    # Infer fluorine position by bond-graph distance from ipso C23 on the phenyl ring C23..C28.
    ring=range(23,29); graph={a:[b for b in ring if b!=a and np.linalg.norm(xyz[a-1]-xyz[b-1])<1.65] for a in ring}
    distances={23:0};queue=[23]
    for a in queue:
        for b in graph[a]:
            if b not in distances:distances[b]=distances[a]+1;queue.append(b)
    assert len(distances)==6 and all(len(v)==2 for v in graph.values())
    fpositions=[]
    for a in f:
        carbon=min(ring,key=lambda b:np.linalg.norm(xyz[a-1]-xyz[b-1]))
        assert np.linalg.norm(xyz[a-1]-xyz[carbon-1])<1.6
        fpositions.append(dict(system=system,F_atom=int(a),bonded_ring_C=int(carbon),position={1:'ortho',2:'meta',3:'para'}[distances[carbon]]))
    # Include all atoms bonded to phenyl C23..C28 outside the core, including H replacements.
    sub=set(ring)
    for i,z in enumerate(zs):
        if z in (1,9) and min(np.linalg.norm(xyz[i]-xyz[c-1]) for c in ring)<1.6:sub.add(i+1)
    amine=list(range(n,n+13));core=list(range(1,23))
    return mapping,dict(amine_side=amine,coumarin_core=core,substituent_side=sorted(sub)),hs,fpositions


def main():
    config=yaml.safe_load((HERE/'config/systems.yaml').read_text());data=ROOT/'data/completed_20261017'
    restore_inputs(data)
    out=HERE/'outputs';csv=out/'csv';fig=out/'figures';logs=out/'oriented_logs';reports=out/'reports'
    for p in (csv,fig,logs,reports):p.mkdir(parents=True,exist_ok=True)
    baseconfig=yaml.safe_load((ROOT/'analysis2_force/config/analysis2_config.yaml').read_text())
    targets=yaml.safe_load((ROOT/'base_analysis/config/analysis_targets.yaml').read_text());targets['fragments']={}
    elements=pd.read_csv(ROOT/'base_analysis/config/atom_map.csv')[['alias','expected_element']]
    atom_map=elements.copy();records={};stationary_rows=[];fpositions=[];h_atoms={};warnings=[]
    lys=parse_stationary(data/'Lys.log');stationary_rows.append(dict(system='Lys',role='reactant',**serializable(lys)))
    for system,spec in config['systems'].items():
        records[system]={role:parse_stationary(data/spec[role]) for role in ('reactant','product','ts')}
        for role,record in records[system].items():
            stationary_rows.append(dict(system=system,role=role,**serializable(record)))
            if not record['normal_termination']:warnings.append(f'{system} {role}: normal termination missing')
            if role!='ts' and record['imaginary_frequencies']:warnings.append(f"{system} {role}: imaginary frequencies {record['imaginary_frequencies']} cm^-1; minimum not fully verified")
            if role=='ts' and len(record['imaginary_frequencies'])!=1:raise ValueError(f'{system}: expected one imaginary TS frequency')
        validate_reference(records[system]['reactant'],records[system]['ts'],records[system]['product'],lys)
        mapping,fragments,hs,fp=mapping_from_ts(system,records[system]['ts'],spec['reaction_N'])
        atom_map[system]=atom_map.alias.map(mapping);targets['fragments'][system]=fragments;fpositions+=fp;h_atoms[system]=hs
        for row in atom_map.itertuples():
            assert records[system]['ts']['atomic_numbers'][mapping[row.alias]-1]=={'C':6,'N':7,'O':8,'F':9,'H':1}[row.expected_element]
    atom_map.to_csv(HERE/'config/atom_map.csv',index=False)
    (HERE/'config/analysis_targets.yaml').write_text(yaml.safe_dump(targets,sort_keys=False))
    table(fpositions,csv/'fluorine_positions.csv');table(stationary_rows,csv/'stationary_point_audit.csv')
    (reports/'stationary_point_audit.json').write_text(json.dumps(stationary_rows,indent=2)+'\n')
    frames={};profiles={};marks={};splines={};orientation=[];energy_rows=[];endpoint_rows=[];geometry_rows=[];dipoles=[]
    validations=[];fragments_all=[];alignments=[];contacts=[];correlations=[];summary_rows=[];sensitivity=[];com_rows=[]
    energies_for_plot={}
    for system,spec in config['systems'].items():
        print(f'{system}: preparing oriented IRC',flush=True)
        stationary=records[system];ts=stationary['ts'];reactant=stationary['reactant'];product=stationary['product'];n=spec['reaction_N']
        frame,audit=prepare_irc(data/spec['irc'],ts,system,n,logs/f'{system}_oriented.log')
        for key in ('source','oriented_log'):audit[key]=str(Path(audit[key]).relative_to(ROOT))
        orientation.append(audit);frames[system]=frame
        if not audit['CN_monotonically_decreasing']:warnings.append(f"{system}: local C-N increases up to {audit['max_CN_increase_A']:.5f} A; overall orientation follows long-to-short endpoints")
        if audit['excluded_unconverged_points']:warnings.append(f"{system}: omitted {len(audit['excluded_unconverged_points'])} unconverged terminal point(s); no extrapolated coordinates used")
        point=frame[['signed_irc_coordinate','path_number','point_number']].drop_duplicates().sort_values('signed_irc_coordinate')
        profile,landmark,spline=reaction_force_profile(frame,KCAL,4001,None)
        profile['E_relative_separated_reactants_kcal_mol']=(profile.E_raw_hartree-reactant['energy_hartree']-lys['energy_hartree'])*KCAL
        profiles[system]=profile;marks[system]=landmark;splines[system]=spline
        profile.to_csv(csv/f'reaction_force_{system}.csv',index=False)
        frame.to_csv(csv/f'parsed_oriented_{system}.csv.gz',index=False,compression={'method':'gzip','mtime':0})
        pd.DataFrame([vars(landmark)]).to_csv(csv/f'landmarks_{system}.csv',index=False)
        for label,lam in [('interpolation',0.0),('GCV',None),('lambda_0.0001',1e-4),('lambda_0.01',.01)]:
            _,lm,sp=reaction_force_profile(frame,KCAL,4001,lam)
            sensitivity.append(dict(system=system,smoothing=label,alpha=lm.alpha_continuous,gamma=lm.gamma_continuous,W2_kcal_mol=(ts['energy_hartree']-float(sp(lm.alpha_continuous)))*KCAL))
        lo=frame[frame.signed_irc_coordinate==frame.signed_irc_coordinate.min()].sort_values('atom_number')
        hi=frame[frame.signed_irc_coordinate==frame.signed_irc_coordinate.max()].sort_values('atom_number')
        _,prod_rmsd=aligned(product['xyz'],hi[['x','y','z']].to_numpy())
        ealpha=float(spline(landmark.alpha_continuous));eleft=float(lo.energy_hartree.iloc[0])
        e=dict(system=system,**energy_terms(reactant,ts,product,lys,ealpha,eleft))
        energy_rows.append(e)
        composition=reactant['composition'].copy()
        for z,count in lys['composition'].items():composition[z]=composition.get(z,0)+count
        assert composition==ts['composition']==product['composition']
        for role,rec in stationary.items():
            assert rec['method_basis']==ts['method_basis'] and rec['solvent']==ts['solvent'],(system,role,'calculation mismatch')
            assert rec['charge']==ts['charge']==0 and rec['multiplicity']==ts['multiplicity']==1
        endpoint_rows.append(dict(system=system,product_opt_vs_IRC_RMSD_A=prod_rmsd,
            product_opt_minus_IRC_endpoint_kcal_mol=(product['energy_hartree']-float(hi.energy_hartree.iloc[0]))*KCAL,
            product_opt_CN_A=float(np.linalg.norm(product['xyz'][19]-product['xyz'][n-1])),IRC_product_CN_A=bond(hi,20,n),
            composition_balanced_with_Lys=True))
        if (product['energy_hartree']-float(hi.energy_hartree.iloc[0]))*KCAL>.5:warnings.append(f'{system}: supplied optimized product is higher than the IRC product-side endpoint; conformer/path correspondence needs checking')
        geo,details=analyze_geometry_context(system,frame,atom_map,targets,baseconfig['geometry_context'])
        geo.to_csv(csv/f'geometry_change_{system}.csv',index=False);geometry_rows.append(dict(system=system,**{k:v for k,v in details.items() if k!='targets'}))
        atomic=analyze_atomic_forces(system,frame,atom_map,targets);atomic.to_csv(csv/f'atomic_force_{system}.csv',index=False)
        ac,validation,fragment,alignment=analyze_force_motion(system,frame,profile,atom_map,targets,baseconfig['force_motion'],KCAL)
        ac.to_csv(csv/f'atomic_contribution_{system}.csv.gz',index=False,compression={'method':'gzip','mtime':0})
        validations.append(validation);fragments_all.append(fragment);alignments.append(alignment)
        plots.energy_force_kappa(system,profile,landmark,spline,reactant['energy_hartree']+lys['energy_hartree'],KCAL,fig/f'01_{system}_energy_force_kappa.png')
        plots.target_geometry(system,frame,atom_map,targets,landmark,fig/f'04_{system}_geometry.png')
        plots.geometry_context(system,geo,landmark,atom_map,targets,fig/f'05_{system}_geometry_context.png')
        plots.atomic_force_context(system,atomic,landmark,atom_map,targets,fig/f'06_{system}_atomic_force.png')
        energies_for_plot[system]={'reactant':SimpleNamespace(energy_hartree=ts['energy_hartree'])}
        mapping=atom_map.set_index('alias')[system].to_dict();system_contacts=[]
        fs=[p['F_atom'] for p in fpositions if p['system']==system]
        for row in point.itertuples(index=False):
            p=frame[(frame.path_number==row.path_number)&(frame.point_number==row.point_number)].sort_values('atom_number');xyz=p[['x','y','z']].to_numpy();r=float(row.signed_irc_coordinate)
            region='reactant' if r<landmark.alpha_continuous else ('transition' if r<=landmark.gamma_continuous else 'product')
            for ff in fs:
                for hh in h_atoms[system]:
                    system_contacts.append(dict(system=system,R=r,region=region,F_atom=ff,H_atom=hh,N_atom=n,
                        F_position=next(t['position'] for t in fpositions if t['system']==system and t['F_atom']==ff),
                        FH_A=float(np.linalg.norm(xyz[ff-1]-xyz[hh-1])),NHF_angle_deg=angle(xyz[n-1],xyz[hh-1],xyz[ff-1]),
                        CN_A=bond(p,20,n),selected_at_TS=ff==mapping['hb_substituent_F'] and hh==mapping['hb_substituent_H']))
            masses=p.atomic_number.map({1:1.008,6:12.011,7:14.007,8:15.999,9:18.998403163}).to_numpy()
            cms={name:np.average(xyz[np.array(ids)-1],axis=0,weights=masses[np.array(ids)-1]) for name,ids in targets['fragments'][system].items()}
            for name,specbase in targets['com_baselines'].items():
                a,b=specbase['fragments'];com_rows.append(dict(system=system,R=r,baseline=name,COM_distance_A=float(np.linalg.norm(cms[a]-cms[b]))))
        c=pd.DataFrame(system_contacts);contacts+=system_contacts
        selected=c[c.selected_at_TS].sort_values('R').copy();selected['dFH_dR']=np.gradient(selected.FH_A,selected.R,edge_order=2)
        selected['force']=profile.F_smooth_kcal_mol_per_R.to_numpy();selected['dCN_dR']=np.gradient(selected.CN_A,selected.R,edge_order=2)
        for region,mask in [('all',np.ones(len(selected),bool)),('reactant',selected.R<landmark.alpha_continuous),('transition',(selected.R>=landmark.alpha_continuous)&(selected.R<=landmark.gamma_continuous)),('pre_TS',selected.R<0),('post_TS',selected.R>0)]:
            sub=selected.loc[mask]
            for a,b in [('FH_A','CN_A'),('FH_A','force'),('dFH_dR','force'),('dFH_dR','dCN_dR')]:
                correlations.append(dict(system=system,region=region,x=a,y=b,n=len(sub),pearson_r=float(pearsonr(sub[a],sub[b]).statistic),spearman_r=float(spearmanr(sub[a],sub[b]).statistic),inferential_p_value_not_used=True))
        tst=selected.loc[selected.R.abs().idxmin()]
        summary_rows.append(dict(system=system,ortho_F_count=sum(p['position']=='ortho' for p in fpositions if p['system']==system),
            TS_CN_A=tst.CN_A,TS_nearest_FH_A=tst.FH_A,TS_NHF_angle_deg=tst.NHF_angle_deg,selected_F_atom=int(tst.F_atom),selected_H_atom=int(tst.H_atom),
            min_selected_FH_A=float(selected.FH_A.min()),reactant_side_FH_A=float(selected.FH_A.iloc[0]),product_side_FH_A=float(selected.FH_A.iloc[-1]),
            alpha=landmark.alpha_continuous,gamma=landmark.gamma_continuous,imaginary_frequency_cm1=ts['imaginary_frequencies'][0]))
        # Original path numbers map to the original dipoles (before orientation export swaps paths).
        dm=parse_irc_dipoles(data/spec['irc'])
        for row in point.itertuples(index=False):
            original_path=(3-row.path_number) if audit['coordinate_multiplier']==-1 and row.path_number else row.path_number
            v=dm.get((int(original_path),int(row.point_number)))
            dipoles.append(dict(system=system,R=row.signed_irc_coordinate,dipole_magnitude_debye=float(np.linalg.norm(v)) if v is not None else np.nan))
        print(f'{system}: {len(profile)} points; C-N {audit["reactant_side_CN_A"]:.3f} -> {audit["product_side_CN_A"]:.3f} A',flush=True)
    table(orientation,csv/'orientation_audit.csv');(reports/'orientation_audit.json').write_text(json.dumps(orientation,indent=2)+'\n')
    energy=pd.DataFrame(energy_rows)
    energy.to_csv(csv/'energy_comparison.csv',index=False)
    summary=pd.DataFrame(summary_rows).merge(energy[['system','absolute_electronic_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol','IRC_left_to_TS_kcal_mol','alpha_to_TS_kcal_mol','RRHO_G_provisional']],on='system');summary.to_csv(csv/'system_summary.csv',index=False)
    # Retire the prior inter-system centering output when regenerating in place.
    (csv/'MoriPF_baseline_differences.csv').unlink(missing_ok=True)
    table(endpoint_rows,csv/'endpoint_correspondence.csv');table(geometry_rows,csv/'geometry_population_summary.csv')
    table(contacts,csv/'all_NH_F_contacts.csv');table(correlations,csv/'descriptive_correlations.csv');table(sensitivity,csv/'smoothing_sensitivity.csv')
    com=table(com_rows,csv/'COM_distances.csv');align=pd.concat(alignments,ignore_index=True)
    corrected=[]
    for (system,target),g in align.groupby(['system','target']):
        baseline=targets['distances'][target]['com_baseline'];base=com[(com.system==system)&(com.baseline==baseline)].set_index('R').COM_distance_A
        tsdist=float(g.loc[g.R.abs().idxmin(),'q_angstrom']);tscom=float(base.loc[0.0])
        for row in g.itertuples():corrected.append(dict(system=system,target=target,R=row.R,raw_distance_A=row.q_angstrom,delta_from_TS_A=row.q_angstrom-tsdist,COM_baseline=baseline,COM_corrected_A=(row.q_angstrom-tsdist)-(base.loc[row.R]-tscom)))
    table(corrected,csv/'COM_corrected_distances.csv')
    validation=pd.concat(validations,ignore_index=True);validation.to_csv(csv/'force_motion_validation.csv',index=False)
    vr=[]
    for system,g in validation.groupby('system',sort=False):
        g=g.dropna();res=g.sum_atomic_C_hartree_per_R-g.reaction_force_hartree_per_R
        vr.append(dict(system=system,n=len(g),pearson_r=float(g[['sum_atomic_C_hartree_per_R','reaction_force_hartree_per_R']].corr().iloc[0,1]),normalized_RMSE=float(np.sqrt(np.mean(res**2)/np.mean(g.reaction_force_hartree_per_R**2)))))
    vsummary=table(vr,csv/'force_motion_validation_summary.csv')
    for row in vr:
        if row['pearson_r']<.95 or row['normalized_RMSE']>.1:warnings.append(f"{row['system']}: force-motion reconstruction deviates substantially; inspect validation before interpreting fragment contributions")
    fragment=pd.concat(fragments_all,ignore_index=True);fragment.to_csv(csv/'fragment_contributions.csv',index=False);align.to_csv(csv/'force_motion_alignment.csv',index=False)
    dipole=table(dipoles,csv/'dipoles.csv')
    for system,stationary in records.items():
        profiles[system].to_csv(csv/f'reaction_force_{system}.csv',index=False)
        offset=stationary['reactant']['energy_hartree']+lys['energy_hartree']
        energies_for_plot[system]={'reactant':SimpleNamespace(energy_hartree=offset)}
    plots.energy_force_kappa_comparison(profiles,marks,splines,energies_for_plot,KCAL,fig/'01_four_systems_energy_force_kappa.png')
    plots.target_geometry_comparison(frames,atom_map,targets,marks,fig/'04_four_systems_geometry.png')
    plots.atomic_force_motion_validation(validation,marks,fig/'08_force_motion_validation.png')
    plots.fragment_contribution_plot(fragment,marks,fig/'09_fragment_contributions.png')
    plots.force_motion_alignment_plot(align,marks,atom_map,targets,fig/'10_force_motion_alignment.png')
    plots.dipole_comparison(dipole,marks,fig/'07_dipoles.png')
    from reporting import create_extra_figures,write_report
    create_extra_figures(summary,energy,pd.DataFrame(contacts),pd.DataFrame(corrected),fig,config)
    write_report(summary,energy,orientation,stationary_rows,pd.DataFrame(endpoint_rows),vsummary,warnings,reports,config)
    from two_sided_contacts import main as analyze_two_sided
    analyze_two_sided()
    from contact_specificity import main as analyze_specificity
    analyze_specificity()
    (reports/'execution.json').write_text(json.dumps(dict(warnings=warnings,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,systems=list(records),energy_reference='each separated substrate + Lys',Lys_sha256=lys['sha256'],Lys_method_basis=lys['method_basis'],Lys_solvent=lys['solvent']),indent=2)+'\n')
    from thermochemistry_profiles import main as analyze_thermochemistry
    analyze_thermochemistry()
    from figure_explanations import write_guide
    write_guide()
    manifest_path=reports/'output_manifest.json'
    manifest=[dict(path=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=sha256(p)) for p in sorted(out.rglob('*')) if p.is_file() and p!=manifest_path]
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    print(summary.to_string(index=False));print('\n'.join(warnings))

if __name__=='__main__':main()
