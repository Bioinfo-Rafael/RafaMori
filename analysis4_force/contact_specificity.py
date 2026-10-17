"""Compare nominated contacts with all atoms and geometry-selected controls.

No electronic-structure calculation. Pair quantities are projections of total
atomic forces, not decompositions of pair interaction energy.
"""
from __future__ import annotations
import itertools,json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from scipy.sparse.csgraph import shortest_path
from scipy.stats import pearsonr,spearmanr
from two_sided_contacts import HERE,KCAL,BOHR_A,_trajectory,_align_to_ts,ortho_substituent

TARGETS={'amine_contact':tuple(sorted(['X24','NH_near'])),
         'opposite_contact':tuple(sorted(['X25','core21']))}
ENERGIES=['absolute_electronic_barrier_kcal_mol','absolute_RRHO_G_barrier_kcal_mol']

def pairkey(a,b):return '|'.join(sorted([a,b]))

def topology(x,z,n,config):
    radii=np.array([config['covalent_radii_A'][int(zz)] for zz in z])
    distance=np.linalg.norm(x[:,None]-x[None,:],axis=2)
    graph=(distance<config['covalent_radius_scale']*(radii[:,None]+radii[None,:]))&(~np.eye(len(z),dtype=bool))
    graph[19,n]=graph[n,19]=False
    # No inferred bond joins the two original molecules.
    amine=set(range(n,n+13))
    for i in range(len(z)):
        for j in range(len(z)):
            if (i in amine)!=(j in amine):graph[i,j]=False
    assert all(graph[i].sum()==1 for i in range(len(z)) if z[i] in [1,9])
    return graph,shortest_path(graph,directed=False,unweighted=True)

def map_sites(system,x,z,n,graph):
    """One-to-one sites only; chemically equivalent alkyl H atoms stay unmapped."""
    labels={i:f'core{i+1:02d}' for i in range(28)}
    for c in range(24,29):labels[ortho_substituent(x,z,c)]=f'X{c}'
    hs=[i for i in np.where(graph[n])[0] if z[i]==1]
    near=next(i for i,v in labels.items() if v=='X24')
    hs=sorted(hs,key=lambda h:np.linalg.norm(x[h]-x[near]))
    assert len(hs)==2
    labels[n]='amineN';labels[hs[0]]='NH_near';labels[hs[1]]='NH_other'
    # Heavy-atom numbering/order shared by the supplied input structures.
    for i,j in enumerate([n-11,n-7,n-4]):
        assert z[j]==6;labels[j]=f'sub_alkylC{i+1}'
    for i,j in enumerate([n+3,n+6,n+9]):
        assert z[j]==6;labels[j]=f'amineC{i+1}'
    assert len(set(labels.values()))==len(labels)
    return labels

def interpolate(x,y,newx):
    if y.ndim==1:return np.interp(newx,x,y)
    return np.stack([np.interp(newx,x,col) for col in y.T],axis=1)

def calculate_system(system,frame,windows,config):
    keys,atoms,xyz,forces=_trajectory(frame);r=keys.signed_irc_coordinate.to_numpy(float)
    ts=int(np.argmin(abs(r)));z=frame[frame.signed_irc_coordinate==r[ts]].sort_values('atom_number').atomic_number.to_numpy(int)
    n=44 if system=='MoriPF' else 42
    graph,gd=topology(xyz[0],z,n,config)
    labels=map_sites(system,xyz[ts],z,n,graph)
    label=[labels.get(i,f'unmapped_{i+1}') for i in range(57)]
    amine=np.array([n<=i<n+13 for i in range(57)])
    edges=[dict(system=system,atom_a=i+1,atom_b=j+1,site_a=label[i],site_b=label[j]) for i,j in zip(*np.where(np.triu(graph,1)))]
    mapping=[dict(system=system,atom=i+1,element={1:'H',6:'C',7:'N',8:'O',9:'F'}[int(z[i])],site=label[i],cross_system_mapped=i in labels,molecule='amine' if amine[i] else 'substrate') for i in range(57)]
    xyz,forces=_align_to_ts(xyz,forces,ts)
    velocity=np.gradient(xyz,r,axis=0,edge_order=2)/BOHR_A
    atomic_C=np.einsum('pij,pij->pi',forces,velocity);norm=np.linalg.norm(forces,axis=2)
    aa,bb=np.triu_indices(57,1)
    delta=xyz[:,bb]-xyz[:,aa];distance=np.linalg.norm(delta,axis=2);unit=delta/distance[:,:,None]
    q=.5*np.einsum('pij,pij->pi',forces[:,bb]-forces[:,aa],unit)
    dr=np.einsum('pij,pij->pi',velocity[:,bb]-velocity[:,aa],unit);power=q*dr
    cn=np.linalg.norm(xyz[:,19]-xyz[:,n],axis=1)
    good=(r<0)&np.isfinite(forces).all(axis=(1,2));rp=r[good];cp=cn[good]
    assert (np.diff(cp)<0).all()
    atoms_out=[];pairs_out=[]
    for window,(start,end) in windows.items():
        assert cp.max()>=start>end>=cp.min()
        r0,r1=np.interp([start,end],cp[::-1],rp[::-1]);rr=np.r_[r0,rp[(rp>r0)&(rp<r1)],r1]
        cc=np.linspace(end,start,501)
        def integ(arr):return np.trapezoid(interpolate(rp,arr[good],rr),rr,axis=0)
        def avg(arr,absolute=False):
            v=interpolate(cp[::-1],arr[good][::-1],cc)
            return np.trapezoid(abs(v) if absolute else v,cc,axis=0)/(start-end)
        aw=-integ(atomic_C)*KCAL;meanforce=avg(norm)*1000
        w=-integ(power)*KCAL;mq=avg(q,True)*1000;md=avg(distance)
        # Exact minimum for the piecewise-linear sampled distance, including boundaries.
        mind=interpolate(rp,distance[good],rr).min(axis=0)
        for i in range(57):
            atoms_out.append(dict(system=system,window=window,**{k:v for k,v in mapping[i].items() if k!='system'},atomic_work_kcal_mol=aw[i],mean_force_mEh_bohr=meanforce[i]))
        for k,(a,b) in enumerate(zip(aa,bb)):
            category='intermolecular' if amine[a]!=amine[b] else ('intramolecular_amine' if amine[a] else 'intramolecular_substrate')
            key=pairkey(label[a],label[b]);nominated=next((name for name,sites in TARGETS.items() if key==pairkey(*sites)),'')
            pairs_out.append(dict(system=system,window=window,pair=key,atom_a=a+1,atom_b=b+1,site_a=label[a],site_b=label[b],element_a=mapping[a]['element'],element_b=mapping[b]['element'],category=category,
                cross_system_mapped=a in labels and b in labels,graph_separation=gd[a,b],min_distance_A=mind[k],mean_distance_A=md[k],
                eligible=bool(gd[a,b]>=config['minimum_graph_separation'] and mind[k]<=config['contact_cutoff_A']),
                nominated=nominated,pair_atomic_work_kcal_mol=aw[a]+aw[b],radial_work_kcal_mol=w[k],mean_abs_projection_mEh_bohr=mq[k]))
    return atoms_out,pairs_out,mapping,edges

def pool_mask(g,target,pool,config):
    mask=g.eligible&(g.category==target.category)
    if pool=='distance_matched':mask &= abs(g.mean_distance_A-target.mean_distance_A)<=config['matched_mean_distance_tolerance_A']
    elif pool=='shared_anchor':
        anchor='core21' if target.nominated=='opposite_contact' else 'NH_near'
        mask &= (g.site_a==anchor)|(g.site_b==anchor)
    elif pool=='disjoint':
        sites={target.site_a,target.site_b}
        mask &= ~g.site_a.isin(sites)&~g.site_b.isin(sites)
    elif pool in ['chemically_related','distance_and_chemistry']:
        if target.nominated=='opposite_contact':
            chemistry=(g.element_a.isin(['H','F'])&g.element_b.isin(['O','N']))|(g.element_b.isin(['H','F'])&g.element_a.isin(['O','N']))
        else:
            chemistry=((g.element_a=='H')&g.element_b.isin(['H','F','O','N']))|((g.element_b=='H')&g.element_a.isin(['H','F','O','N']))
        mask &= chemistry
        if pool=='distance_and_chemistry':mask &= abs(g.mean_distance_A-target.mean_distance_A)<=config['matched_mean_distance_tolerance_A']
    elif pool=='exclude_1_4':mask &= g.graph_separation>=4
    elif pool!='same_category':raise ValueError(pool)
    return mask

POOLS=['same_category','distance_matched','shared_anchor','disjoint','chemically_related','distance_and_chemistry','exclude_1_4']

def within_rankings(atoms,pairs,config):
    rows=[]
    for (system,window),g in atoms.groupby(['system','window'],sort=False):
        for site in ['X24','NH_near','X25','core21']:
            t=g[g.site==site].iloc[0]
            for metric in config['atom_metrics']:
                v=abs(t[metric]);vals=abs(g[metric])
                rows.append(dict(kind='atom',system=system,window=window,target=site,pool='all_57_atoms',metric=metric,value=t[metric],magnitude=v,n_controls=len(g)-1,stronger_controls=int((vals>v+1e-12).sum()),magnitude_rank=int((vals>v+1e-12).sum())+1,median_control_magnitude=float(vals[g.site!=site].median())))
    for (system,window),g in pairs.groupby(['system','window'],sort=False):
        for _,t in g[g.nominated!=''].iterrows():
            assert t.eligible
            for pool in POOLS:
                controls=g[pool_mask(g,t,pool,config)&(g.pair!=t.pair)]
                for metric in config['pair_metrics']:
                    v=abs(t[metric]);vals=abs(controls[metric])
                    rows.append(dict(kind='contact',system=system,window=window,target=t.nominated,pool=pool,metric=metric,value=t[metric],magnitude=v,n_controls=len(controls),stronger_controls=int((vals>v+1e-12).sum()),magnitude_rank=int((vals>v+1e-12).sum())+1,median_control_magnitude=float(vals.median())))
    return pd.DataFrame(rows)

def correlation_rows(data,key,metrics,energy,kind):
    rows=[]
    systems=energy.index.tolist()
    for (window,label),g in data.groupby(['window',key],sort=False):
        if len(g)!=4 or set(g.system)!=set(systems):continue
        g=g.set_index('system').loc[systems]
        for metric,target in itertools.product(metrics,ENERGIES):
            x=g[metric].to_numpy();y=energy[target].to_numpy()
            if np.ptp(x)<1e-12:continue
            r=float(pearsonr(x,y).statistic);rho=float(spearmanr(x,y).statistic)
            leave=[float(pearsonr(np.delete(x,i),np.delete(y,i)).statistic) for i in range(4)]
            rows.append(dict(kind=kind,window=window,item=label,metric=metric,barrier=target,pearson_r=r,spearman_r=rho,leave_one_out_r_min=min(leave),leave_one_out_r_max=max(leave),median_magnitude=float(np.median(abs(x)))))
    return pd.DataFrame(rows)

def specificity(pairs,corr,energy,config):
    summaries=[];examples=[];permrows=[]
    systems=energy.index.tolist()
    eligible=pairs[pairs.eligible & pairs.cross_system_mapped]
    for window,g in eligible.groupby('window',sort=False):
        common=set(g.groupby('pair').filter(lambda f:len(f)==4).pair)
        g=g[g.pair.isin(common)]
        for targetname,sites in TARGETS.items():
            key=pairkey(*sites);t=g[g.pair==key].set_index('system').loc[systems]
            for pool in POOLS:
                controlsets=[]
                for system in systems:
                    local=g[g.system==system];row=t.loc[system]
                    controlsets.append(set(local.loc[pool_mask(local,row,pool,config),'pair'])-{key})
                keys=set.intersection(*controlsets)
                for metric,barrier in itertools.product(config['pair_metrics'],ENERGIES):
                    c=corr[(corr.kind=='contact')&(corr.window==window)&(corr.metric==metric)&(corr.barrier==barrier)].set_index('item')
                    keys_valid=sorted(keys&set(c.index));controls=c.loc[keys_valid];tr=c.loc[key]
                    tv=t[metric].to_numpy();threshold=abs(tr.pearson_r)
                    simultaneous=[]
                    for item in keys_valid:
                        vals=g[g.pair==item].set_index('system').loc[systems,metric].to_numpy()
                        if (abs(vals)>=abs(tv)-1e-12).all():simultaneous.append(item)
                    strong=controls.loc[simultaneous]
                    summaries.append(dict(window=window,target=targetname,pool=pool,metric=metric,barrier=barrier,
                        target_pearson_r=tr.pearson_r,target_spearman_r=tr.spearman_r,n_common_controls=len(controls),
                        controls_at_least_as_correlated=int((abs(controls.pearson_r)>=threshold-1e-12).sum()),
                        abs_correlation_rank=int((abs(controls.pearson_r)>threshold+1e-12).sum())+1,
                        stronger_in_all_four=len(strong),stronger_and_at_least_as_correlated=int((abs(strong.pearson_r)>=threshold-1e-12).sum())))
                    # Retain strongest competing correlations, including zero/both shared atoms.
                    ordered=controls.assign(absr=abs(controls.pearson_r)).sort_values('absr',ascending=False)
                    dominant=ordered.loc[ordered.index.isin(simultaneous)&(ordered.absr>=threshold-1e-12)].head(3)
                    competitors=pd.concat([ordered.head(5),dominant]).loc[lambda df:~df.index.duplicated()]
                    for item,row in competitors.iterrows():
                        vals=g[g.pair==item].set_index('system').loc[systems]
                        examples.append(dict(window=window,target=targetname,pool=pool,metric=metric,barrier=barrier,control_pair=item,pearson_r=row.pearson_r,spearman_r=row.spearman_r,stronger_in_all_four=item in simultaneous,
                            **{s+'_value':vals.loc[s,metric] for s in systems},
                            **{s+'_atoms':f'{vals.loc[s,"atom_a"]}-{vals.loc[s,"atom_b"]}' for s in systems}))
                    # Exhaust all 24 label permutations, simultaneously over every control.
                    # Descriptive scan calibration, NOT a confirmatory significance test.
                    matrix=g[g.pair.isin(keys_valid+[key])].pivot(index='pair',columns='system',values=metric).reindex(columns=systems).to_numpy(copy=True)
                    matrix-=matrix.mean(axis=1,keepdims=True);matrix/=np.linalg.norm(matrix,axis=1,keepdims=True)
                    y=energy[barrier].to_numpy();y=y-y.mean();y=y/np.linalg.norm(y)
                    maxima=[float(np.max(abs(matrix@yp))) for yp in itertools.permutations(y)]
                    permrows.append(dict(window=window,target=targetname,pool=pool,metric=metric,barrier=barrier,n_label_permutations=24,observed_target_abs_r=threshold,
                        permutations_with_scan_max_at_least_target=int(sum(v>=threshold-1e-12 for v in maxima)),confirmatory_p_value=False))
    return pd.DataFrame(summaries),pd.DataFrame(examples),pd.DataFrame(permrows)

def main():
    config=yaml.safe_load((HERE/'config/contact_specificity.yaml').read_text())
    windows=yaml.safe_load((HERE/'config/two_sided_contacts.yaml').read_text())
    csv=HERE/'outputs/csv';energy=pd.read_csv(csv/'energy_comparison.csv').set_index('system')
    aa=[];pp=[];mm=[];ee=[]
    for system in energy.index:
        a,p,m,e=calculate_system(system,pd.read_csv(csv/f'parsed_oriented_{system}.csv.gz'),windows['windows'],config)
        aa+=a;pp+=p;mm+=m;ee+=e
    atoms=pd.DataFrame(aa);pairs=pd.DataFrame(pp);mapping=pd.DataFrame(mm);edges=pd.DataFrame(ee)
    # Verify equivalent-site heavy-atom/topological mapping before correlating.
    mapped=edges[~edges.site_a.str.startswith('unmapped_')&~edges.site_b.str.startswith('unmapped_')]
    edge_sets=[{pairkey(a,b) for a,b in zip(g.site_a,g.site_b)} for _,g in mapped.groupby('system')]
    assert all(x==edge_sets[0] for x in edge_sets)
    ranks=within_rankings(atoms,pairs,config)
    ac=correlation_rows(atoms[atoms.cross_system_mapped],'site',config['atom_metrics'],energy,'atom')
    pc=correlation_rows(pairs[pairs.cross_system_mapped & pairs.eligible],'pair',config['pair_metrics'],energy,'contact')
    corr=pd.concat([ac,pc],ignore_index=True)
    specificity_table,examples,permutations=specificity(pairs,corr,energy,config)
    for name,df in [('atoms',atoms),('pairs',pairs),('mapping',mapping),('covalent_edges',edges),('within_ranks',ranks),('correlations',corr),('specificity',specificity_table),('counterexamples',examples),('label_permutations',permutations)]:
        df.to_csv(csv/f'contact_controls_{name}.csv',index=False)
    from contact_specificity_report import create_report
    create_report(atoms,pairs,ranks,corr,specificity_table,examples,permutations,config,windows,energy)
    print(specificity_table[(specificity_table.window==windows['primary_window'])&(specificity_table.metric=='radial_work_kcal_mol')&(specificity_table.barrier==ENERGIES[0])].to_string(index=False))

if __name__=='__main__':main()
