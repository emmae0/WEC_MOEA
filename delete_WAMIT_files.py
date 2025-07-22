####delete old WAMIT files
import os

def delete_WAMIT_files(path):

  os.chdir(path)
  
  if os.path.isfile('mpswec.1'):
    os.remove('mpswec.1')
  if os.path.isfile('mpswec.2'):
    os.remove('mpswec.2')
  if os.path.isfile('mpswec.3'):
    os.remove('mpswec.3')
  if os.path.isfile('mpswec.4'):
    os.remove('mpswec.4')
  if os.path.isfile('mpswec.frc'):
    os.remove('mpswec.frc')
  if os.path.isfile('mpswec.pot'):
    os.remove('mpswec.pot')
  if os.path.isfile('mpswec.gdf'):
    os.remove('mpswec.gdf')
  if os.path.isfile('mpswec.out'):
    os.remove('mpswec.out')
  if os.path.isfile('mpswec.p2f'):
    os.remove('mpswec.p2f')
  if os.path.isfile('fnames.wam'):
    os.remove('fnames.wam')
  if os.path.isfile('errorf.log'):
    os.remove('errorf.log')
  if os.path.isfile('errorp.log'):
    os.remove('errorp.log')

