
'''
    Tumor Landscape Analysis (TLA):
    ##############################

        This script reads lines from a study set table
        Each line has parameters of a particular study

        This formating helps faster and more consistent running of the TLA


'''

# %% Imports
import os
import psutil
import sys
import time
import pandas as pd
from argparse import ArgumentParser

from tla_functions import filexists

__version__  = "2.0.2"


# %% Private classes

class Study:
    
  def __init__(self, study, main_pth):
      
      # loads arguments for this study
      self.name = study['name']
      self.raw_path = os.path.join(main_pth, study['raw_path'])
      
      samples_file = filexists(os.path.join(self.raw_path,
                                            study['raw_samples_table']))
      self.samples = pd.read_csv(samples_file)
      self.samples.fillna('', inplace=True)

      self.dat_pth = os.path.join(main_pth, study['data_path'])
      if not os.path.exists(self.dat_pth):
          print("ERROR: data folder " + self.dat_pth + " does not exist!")
          print("... Please run << TLA_region >> for this study!")
          sys.exit()
          
      # reduces classes df to just the accepted types (i.e. `drop=False`)
      classes_file = filexists(os.path.join(self.raw_path,
                                            study['raw_reg_classes_table']))
      classes = pd.read_csv(classes_file)
      classes.drop(classes.loc[classes['drop']].index, inplace=True)
      classes.drop(columns=['drop'], inplace=True)
      classes.reset_index(inplace=True, drop=True)
      classes['class'] = classes['class'].astype(str)
      self.classes = classes
      
      # creates samples table for output
      self.samples_out = pd.DataFrame()
      self.reg_stats = pd.DataFrame()
      self.cls_stats = pd.DataFrame()


  def mergeTables(self):
      
      import glob
      
      dat_pth = self.dat_pth
      
      # loops over slides in this study
      for i, slide in self.samples.iterrows():
          
          # find samples associated to each slide
          # based on existing folders created by 'tla_region_sample.py'
          slideid = slide.sample_ID
          pth = os.path.join(dat_pth, 'results', 'samples')
          aux = glob.glob(os.path.join(pth, slideid + '*/'))
          sids = [f.split('/')[-2] for f in aux]
          del aux
          
          if (len(sids) > 0):
              for sid in sids:
                  # read table for this sample
                  res_pth = os.path.join(pth, sid)
                  reg_pth = os.path.join(pth, sid, 'regs')
                  
                  fsamples = os.path.join(res_pth, sid + '_samples_reg.csv')
                  fstats = os.path.join(res_pth, sid + '_reg_tbl.csv')
                  fregs = os.path.join(reg_pth, 
                                       sid + '_reg_landscape_metrics.csv')
                  fclss = os.path.join(reg_pth, 
                                       sid + '_reg_class_metrics.csv')
                  
                  isvalid = os.path.exists(fsamples) and \
                      os.path.exists(fstats) and \
                          os.path.exists(fregs)
                          
                  if isvalid:
                      # concat sample data to study samples table
                      aux = pd.read_csv(fsamples)
                      self.samples_out = pd.concat([self.samples_out,
                                                    aux],
                                                   ignore_index=True)
                      
                      # concat regional stats
                      aux = pd.read_csv(fregs)
                      aux['sample_ID'] = sid
                      auy = aux.pivot(index='sample_ID', 
                                      columns='metric', 
                                      values='value').reset_index()
                      aux = pd.read_csv(fstats)
                      aux = aux[['sample_ID', 
                                 'adjacency_index', 
                                 'Simpson']]
                      aux = aux.rename(columns={'Simpson':
                                                'simpson_diversity_index'})
                      auy = pd.merge(auy, aux, how='left',
                                     on=['sample_ID']).fillna(0)
                      self.reg_stats = pd.concat([self.reg_stats,
                                                  auy], 
                                                 ignore_index=True)
                      
                      # concat class-level stats 
                      aux = pd.read_csv(fclss)
                      aux.insert(0, 'sample_ID', sid)
                      self.cls_stats = pd.concat([self.cls_stats,
                                                  aux],
                                                 ignore_index=True)
                      
                  else:
                      print("WARNING: " + sid + " dropped from TLA!")
          else:
              print("WARNING: " + slideid + " dropped from TLA!")
                
      # saves study tables
      self.samples_out.sort_values(by=['sample_ID'])
      self.samples_out['index'] = range(0, len(self.samples_out))
      
      f = os.path.join(dat_pth, self.name + '_samples_reg.csv')
      self.samples_out.to_csv(f, index=False, header=True)
          
      f = os.path.join(dat_pth, 'results', 
                       self.name + '_reg_landscape_stats.csv')
      self.reg_stats.to_csv(f, index=False, header=True)
      
      f = os.path.join(dat_pth, 'results', 
                       self.name + '_reg_classes_stats.csv')
      self.cls_stats.to_csv(f, index=False, header=True)


# %% Private Functions

def memuse():
    m = psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2
    return(round(m, 2))


    # %% Main function

def main(args):
    """
    *******  Main function  *******

    """

    # %% start, checks how the program was launched 
    debug = False
    try:
        args
    except NameError:
        #  if not running from the CLI, run in debug mode
        debug = True
    
    # tracks time and memory 
    start = time.time()
    
    if debug:
        # running from the IDE
        # path of directory containing this script
        main_pth = os.path.dirname(os.getcwd())
        argsfile = os.path.join(main_pth, 'pathAI.csv')
    else:
        # running from the CLI using the bash script
        # path to working directory (above /scripts)
        main_pth = os.getcwd()
        argsfile = os.path.join(main_pth, args.argsfile) 

    print("==> The working directory is: " + main_pth)
    
    if not os.path.exists(argsfile):
        print("ERROR: The specified argument file does not exist!")
        sys.exit()
        
    # only the first study in the argument table will be used
    study = Study( pd.read_csv(argsfile).iloc[0], main_pth)
    
    # %% Summarize stats in study tables
    study.mergeTables()
    
    trun = time.strftime('%H:%M:%S', time.gmtime(time.time()-start))
    print('==> TLA-Sum finished. Time elapsed: ', trun, '[HH:MM:SS]')
    print("==> Max memory used: " + str(memuse()) + "[MB]")
    
    #%%
    return(0)        
        

# %% Argument parser
if __name__ == "__main__":

    # Create the parser
    my_parser = ArgumentParser(prog="tla_sum",
                               description="# Sum Processing module for " +
                               "Tumor Landscape Analysis #",
                               allow_abbrev=False)

    # Add the arguments
    my_parser.add_argument('-v', '--version', 
                           action='version', 
                           version="%(prog)s " + __version__)
    
    my_parser.add_argument('argsfile',
                           metavar="argsfile",
                           type=str,
                           help="Argument table file (.csv) for study set")

    # passes arguments object to main
    main(my_parser.parse_args())
