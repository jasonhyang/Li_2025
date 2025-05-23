import pandas as pd
import multiprocessing as mp
import os

import cobra
from cobra.io import load_matlab_model, save_matlab_model
from cobra.sampling.sampling import sample


def itercall_pysampling(samp_pars):
    """
    Function called by iterate_pysampling in order to simulatenously sample multiple
    genome-scale metabolic models at once using multiprocessing.

    """
    
    imat_model = load_matlab_model(f'{samp_pars[1]}/{samp_pars[0]}/{samp_pars[0]}.mat')
    s = sample(model=imat_model, n=samp_pars[2], processes=1)

    s_len = len(s)
    s.loc['Mean'] = s.iloc[0:s_len].mean()
    s.loc['Median'] = s.iloc[0:s_len].median()
    s.loc['SD'] = s.iloc[0:s_len].std()
    s.loc['1st Q'] = s.iloc[0:s_len].quantile(q=0.25)
    s.loc['3rd Q'] = s.iloc[0:s_len].quantile(q=0.75)
    s.loc['IQR'] = s.iloc[0:s_len].quantile(q=0.75) - s.iloc[0:s_len].quantile(q=0.25)
    s.loc['Min'] = s.iloc[0:s_len].min()
    s.loc['Max'] = s.iloc[0:s_len].max()
    s = s.transpose()
    s.to_csv(f'{samp_pars[1]}/{samp_pars[0]}/{samp_pars[0]} Sampling Results.tsv', sep='\t')

    stats = s.iloc[:, -8:]
    stats.columns = pd.MultiIndex.from_tuples(
        [(samp_pars[0], 'Mean'), (samp_pars[0], 'Median'), (samp_pars[0], 'SD'),
         (samp_pars[0], '1st Q'), (samp_pars[0], '3rd Q'), (samp_pars[0], 'IQR'), (samp_pars[0], 'Min'),
         (samp_pars[0], 'Max')])
    stats = stats.reset_index(drop=True)

    return [samp_pars[0], stats]

def iterate_pysampling(gem_dir, exp_dir, num_samp, order_list, condition):
    """
    This function samples multiple genome-scale metabolic models that have been
    outputted by iterate_imat.R in a concurrent manner and outputs two tab-delimited
    files for each model containing the flux sampling results and the summary statistics
    for each reaction.

    Notes
    -----
    Currently, the function is set to automatically use Gurobi as the default solver
    and optGpSampler as the default sampling algorithm. The function automatically
    uses multiprocessing to simultaneously sample each model.

    Parameters
    ----------
    gem_dir : str
        The string directory of the relevant generic genome-scale metabolic model.
    exp_dir : str
        The string directory where iterate_imat.R exported the genome-scale metabolic
        models that need to be sampled.
    num_samp : int
        The integer number of samples that optGpSampler should produce.
    order_list : list
        Ordered list of RNA-Seq sample names automatically produced by iterate_imat.R
        used to organize resulting sampling results.
    condition : str
        The string descriptor to be used to name the resulting sampling result files.

    Outputs
    -------
    'condition_modelname Sampling Results.tsv'
        The tab-delimited file containing the flux sampling results. Last 8 columns
        include summary statistics.
    'condition_modelname SS.tsv'
        The tab-delimited file containing only the summary statistics of the flux
        sampling results along with additional reaction information including
        reaction ID, name, formula, and subsystem.

    """

    cobra_config = cobra.Configuration()
    cobra_config.solver = 'gurobi'
    gem_model = load_matlab_model(gem_dir)

    pars_list = [[samp, exp_dir, num_samp] for samp in os.listdir(exp_dir)]

    with mp.Pool(processes=mp.cpu_count()-1) as p:
        stats_list = p.map(itercall_pysampling, pars_list)

    samp_stats = [stats_df[1] for stats_df in sorted(stats_list, key=lambda x: order_list.index(x[0]))]

    id_list = [reaction.id for reaction in gem_model.reactions]
    id_df = pd.DataFrame(id_list)
    id_df.columns = ['Reaction ID']
    id_df.columns = pd.MultiIndex.from_tuples([('', 'Reaction ID')])

    name_list = [reaction.name for reaction in gem_model.reactions]
    name_df = pd.DataFrame(name_list)
    name_df.columns = ['Reaction Name']
    name_df.columns = pd.MultiIndex.from_tuples([('', 'Reaction Name')])

    rxn_list = [reaction.reaction for reaction in gem_model.reactions]
    rxn_df = pd.DataFrame(rxn_list)
    rxn_df.columns = ['Reaction Formula']
    rxn_df.columns = pd.MultiIndex.from_tuples([('', 'Reaction Formula')])

    sub_list = [reaction.subsystem for reaction in gem_model.reactions]
    sub_df = pd.DataFrame(sub_list)
    sub_df.columns = ['Subsystem']
    sub_df.columns = pd.MultiIndex.from_tuples([('', 'Subsystem')])

    sampsummary_list = [id_df, name_df, rxn_df, sub_df] + samp_stats
    sampsummary_df = pd.concat(sampsummary_list, axis=1)
    sampsummary_df.to_csv(f'{exp_dir}/{condition} SS.tsv', sep='\t', index=False)

def individual_pysampling(exp_dir, num_samp, condition):
    """
    This function samples a single genome-scale metabolic model outputted by
    iterate_imat.R and outputs two tab-delimited files containing the flux sampling
    results and the summary statistics for each reaction.

    Notes
    -----
    Currently, the function is set to automatically use Gurobi as the default solver
    and optGpSampler as the default sampling algorithm. The function automatically
    uses multiple processors for better performance.

    Parameters
    ----------
    exp_dir : str
        The string directory where iterate_imat.R exported the genome-scale metabolic
        model that needs to be sampled.
    num_samp : int
        The integer number of samples that optGpSampler should produce.
    condition : str
        The string descriptor to be used to name the resulting sampling result files.

    Outputs
    -------
    'condition_modelname Sampling Results.tsv'
        The tab-delimited file containing the flux sampling results. Last 8 columns
        include summary statistics.
    'condition_modelname SS.tsv'
        The tab-delimited file containing only the summary statistics of the flux
        sampling results along with additional reaction information including
        reaction ID, name, formula, and subsystem.

    """

    cobra_config = cobra.Configuration()
    cobra_config.solver = 'gurobi'

    samp_name = [i for i in os.listdir(exp_dir) if os.path.isdir(os.path.join(exp_dir, i))][0]

    imat_model = load_matlab_model(f'{exp_dir}/{samp_name}/{samp_name}.mat')
    s = sample(model=imat_model, n=num_samp, processes=os.cpu_count()-1)

    s_len = len(s)
    s.loc['Mean'] = s.iloc[0:s_len].mean()
    s.loc['Median'] = s.iloc[0:s_len].median()
    s.loc['SD'] = s.iloc[0:s_len].std()
    s.loc['1st Q'] = s.iloc[0:s_len].quantile(q=0.25)
    s.loc['3rd Q'] = s.iloc[0:s_len].quantile(q=0.75)
    s.loc['IQR'] = s.iloc[0:s_len].quantile(q=0.75) - s.iloc[0:s_len].quantile(q=0.25)
    s.loc['Min'] = s.iloc[0:s_len].min()
    s.loc['Max'] = s.iloc[0:s_len].max()
    s = s.transpose()
    s.to_csv(f'{exp_dir}/{samp_name}/{samp_name} Sampling Results.tsv', sep='\t')

    stats = s.iloc[:, -8:]
    stats.columns = pd.MultiIndex.from_tuples(
        [(samp_name, 'Mean'), (samp_name, 'Median'), (samp_name, 'SD'),
         (samp_name, '1st Q'), (samp_name, '3rd Q'), (samp_name, 'IQR'), (samp_name, 'Min'),
         (samp_name, 'Max')])
    stats = stats.reset_index(drop=True)

    id_list = [reaction.id for reaction in imat_model.reactions]
    id_df = pd.DataFrame(id_list)
    id_df.columns = ['Reaction ID']
    id_df.columns = pd.MultiIndex.from_tuples([('', 'Reaction ID')])

    name_list = [reaction.name for reaction in imat_model.reactions]
    name_df = pd.DataFrame(name_list)
    name_df.columns = ['Reaction Name']
    name_df.columns = pd.MultiIndex.from_tuples([('', 'Reaction Name')])

    rxn_list = [reaction.reaction for reaction in imat_model.reactions]
    rxn_df = pd.DataFrame(rxn_list)
    rxn_df.columns = ['Reaction Formula']
    rxn_df.columns = pd.MultiIndex.from_tuples([('', 'Reaction Formula')])

    sub_list = [reaction.subsystem for reaction in imat_model.reactions]
    sub_df = pd.DataFrame(sub_list)
    sub_df.columns = ['Subsystem']
    sub_df.columns = pd.MultiIndex.from_tuples([('', 'Subsystem')])

    sampsummary_list = [id_df, name_df, rxn_df, sub_df, stats]
    sampsummary_df = pd.concat(sampsummary_list, axis=1)
    sampsummary_df.to_csv(f'{exp_dir}/{condition}_{samp_name} SS.tsv', sep='\t', index=False)

def summary_sampling(gem_model, num_samp, condition):
    """
    This function samples an already-defined COBRApy genome-scale metabolic model and
    outputs two tab-delimited files containing the flux sampling results and the
    summary statistics for each reaction along with a matlab file of the model.

    Notes
    -----
    Currently, the function is set to automatically use Gurobi as the default solver
    and optGpSampler as the default sampling algorithm. The function automatically
    uses multiple processors for better performance.

    Parameters
    ----------
    gem_model : cobra.Model
        The model from which to sample flux distributions.
    num_samp : int
        The integer number of samples that optGpSampler should produce.
    condition : str
        The string descriptor to be used to name the resulting sampling result files.

    Outputs
    -------
    'condition_modelname.mat'
        The matlab file output for the genome-scale metabolic model.
    'condition_modelname Sampling Results.tsv'
        The tab-delimited file containing the flux sampling results. Last 8 columns
        include summary statistics.
    'condition_modelname SS.tsv'
        The tab-delimited file containing only the summary statistics of the flux
        sampling results along with additional reaction information including
        reaction ID, name, formula, and subsystem.

    """

    cobra_config = cobra.Configuration()
    cobra_config.solver = 'gurobi'
    samp_name = gem_model.id
    save_matlab_model(gem_model, f'{condition}_{samp_name}.mat')

    s = sample(model=gem_model, n=num_samp, processes=os.cpu_count()-1)
    s_len = len(s)
    s.loc['Mean'] = s.iloc[0:s_len].mean()
    s.loc['Median'] = s.iloc[0:s_len].median()
    s.loc['SD'] = s.iloc[0:s_len].std()
    s.loc['1st Q'] = s.iloc[0:s_len].quantile(q=0.25)
    s.loc['3rd Q'] = s.iloc[0:s_len].quantile(q=0.75)
    s.loc['IQR'] = s.iloc[0:s_len].quantile(q=0.75) - s.iloc[0:s_len].quantile(q=0.25)
    s.loc['Min'] = s.iloc[0:s_len].min()
    s.loc['Max'] = s.iloc[0:s_len].max()
    s = s.transpose()
    s.to_csv(f'{condition}_{samp_name} Sampling Results.tsv', sep='\t')

    stats = s.iloc[:, -8:]
    stats.columns = pd.MultiIndex.from_tuples(
        [(samp_name, 'Mean'), (samp_name, 'Median'), (samp_name, 'SD'),
         (samp_name, '1st Q'), (samp_name, '3rd Q'), (samp_name, 'IQR'), (samp_name, 'Min'),
         (samp_name, 'Max')])
    stats = stats.reset_index(drop=True)

    id_list = [reaction.id for reaction in gem_model.reactions]
    id_df = pd.DataFrame(id_list)
    id_df.columns = ['Reaction ID']
    id_df.columns = pd.MultiIndex.from_tuples([('', 'Reaction ID')])

    name_list = [reaction.name for reaction in gem_model.reactions]
    name_df = pd.DataFrame(name_list)
    name_df.columns = ['Reaction Name']
    name_df.columns = pd.MultiIndex.from_tuples([('', 'Reaction Name')])

    rxn_list = [reaction.reaction for reaction in gem_model.reactions]
    rxn_df = pd.DataFrame(rxn_list)
    rxn_df.columns = ['Reaction Formula']
    rxn_df.columns = pd.MultiIndex.from_tuples([('', 'Reaction Formula')])

    sub_list = [reaction.subsystem for reaction in gem_model.reactions]
    sub_df = pd.DataFrame(sub_list)
    sub_df.columns = ['Subsystem']
    sub_df.columns = pd.MultiIndex.from_tuples([('', 'Subsystem')])

    sampsummary_list = [id_df, name_df, rxn_df, sub_df, stats]
    sampsummary_df = pd.concat(sampsummary_list, axis=1)
    sampsummary_df.to_csv(f'{condition}_{samp_name} SS.tsv', sep='\t', index=False)
