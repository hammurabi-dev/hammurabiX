# -*- coding: utf-8 -*-
'''
@brief: python wrapper for c++ routines
@author: Jiaxin Wang
@email: jiwang@sissa.it

@methods:
# Import class
In []: import wrapper as wp

# Initialize object
# exe_name: the name of executable, 'hammurabi' by default
In []: object = wp.Wrapper(exe_name='')

# Modify parameter value from base xml file to temp xml file
In []: object.mod_par(keys=['key1','key2',...],tag='',attrib='')

# Run the executable, read in output maps, delete all temporary files
In []: object.call()
'''

import os
import subprocess
import healpy as hp
import xml.etree.ElementTree as et
import random as rnd
import numpy as np
import tempfile as tf

class Wrapper(object):
	
	# we try to keep __init__ as light as possible
	def __init__(self,exe_name=None):
		# we offer a freedom for different executables
		self.target = 'hammurabi'
		if exe_name is not None:
			self.target = exe_name
		
		# for simplicity we nail down the working directory
		self.dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
		self.executable = os.path.join(self.dir,self.target)
		# nail down parameter file           
		self.base_file = os.path.join(self.dir,'params.xml')
		# re-assigned by new_params_copy
		self.temp_file = self.base_file
		# simulation output
		self.sim_map_name = {}
		self.sim_map = {}
		
	
	def call(self):
		import sys
		import time
		# create new temp parameter file
		if self.temp_file is self.base_file:
			self._new_params_copy()
		
		temp_process = subprocess.Popen([self.executable,self.temp_file],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
		temp_process.wait()
        	'''
        	#while temp_process.poll() is None:
            		#sys.stdout.write('.')
            		#sys.stdout.flush()
    			#time.sleep(1.0)
        	'''
		if temp_process.returncode != 0:
			last_call_log,last_call_err = temp_process.communicate()
			print last_call_log
			print last_call_err
		
		self._get_sims()
		self._del_params_copy()
		
    
	def mod_par(self,keys=None,tag=None,attrib=None):
		if (keys is None or tag is None or attrib is None):
			print 'ERR: missing input'
			exit(1)
		
		# create new temp parameter file
		if self.temp_file is self.base_file:
			self._new_params_copy()
		
		tree = et.parse(self.temp_file)
		root = tree.getroot()
		path_str = '.'
		for key in keys:
			path_str += '/' + key
    		
		target = root.find(path_str)
		if target is None:
			print 'ERR: wrong element path'
			exit(1)
		
		target.set(tag,attrib)
		tree.write(self.temp_file)
    
	
	
	def _get_sims(self):
        	# locate the files
		iqu_path = os.path.join(self.dir,self.sim_map_name['iqu'])
		fd_path = os.path.join(self.dir,self.sim_map_name['fd'])
		dm_path = os.path.join(self.dir,self.sim_map_name['dm'])
		
		if(os.path.isfile(dm_path)):
			[DM] = self._read_fits_file(dm_path)
			self.sim_map['DM'] = DM
			os.remove(dm_path)
		
		if(os.path.isfile(iqu_path)):
			[Is,Qs,Us] = self._read_fits_file(iqu_path)
            		#self.sim_map['I'] = Is
			self.sim_map['Q'] = Qs
			self.sim_map['U'] = Us
			os.remove(iqu_path)
		
		if(os.path.isfile(fd_path)):
			[Fd] = self._read_fits_file(fd_path)
			self.sim_map['F'] = Fd
          		os.remove(fd_path)
	
	
	def _read_fits_file(self,path):
		rslt = []
		i = 0
		while True:
			try:
				loaded_map = hp.read_map(path,verbose=False,field=i)
				rslt += [loaded_map]
				i += 1
			
			except IndexError:
				break
		
           	return rslt



	def _new_params_copy(self):
		# create a random file name which doesn't exist currently
		fd,new_path = tf.mkstemp(prefix='params_',suffix='.xml',dir=self.dir)
		os.close(fd)
		rnd_idx = new_path[new_path.index('params_')+7:-4]
		self.temp_file = new_path
		# copy base_file to temp_file
		tree = et.parse(self.base_file)
		root = tree.getroot()
		self.sim_map_name['iqu'] = 'iqu_sync_'+rnd_idx+'.fits'
		self.sim_map_name['fd'] = 'fd_'+rnd_idx+'.fits'
		self.sim_map_name['dm'] = 'dm_'+rnd_idx+'.fits'
		root.find('./Output/Sync').set('filename',self.sim_map_name['iqu'])
		root.find('./Output/Faraday').set('filename',self.sim_map_name['fd'])
		root.find('./Output/DM').set('filename',self.sim_map_name['dm'])
		# automatically create a new file
		tree.write(self.temp_file)
	

	
	def _del_params_copy(self):
		if self.temp_file is self.base_file:
			print 'ERR: _del_params_copy failed'
			exit(1)
		else:
			os.remove(self.temp_file)
			self.temp_file = self.base_file
			
		
		
