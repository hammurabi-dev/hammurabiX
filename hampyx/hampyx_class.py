# -*- coding: utf-8 -*-
'''
Developed by Joe Taylor, based on the 
initial work by Theo Steininger

temporary/testing version of hammurabiX python wrapper

warning:

this wrapper file must be copied to hammurabiX executable directory
also it will pass data through disk

methods:

# Import class
In []: import hampyx_class as hamx

# Initialize object
# exe_name: the name of executable, 'hamx' by default
In []: object = hamx.Hampy(exe_type='')

# Modify parameter value from base xml file to temp xml file
In []: object.mod_par(keys=['key1','key2',...],tag='',attrib='')
The strings 'key1', 'key2', etc represent the path to the desired
parameter, going through the xml.

The "tag" is the label for the parameter: eg. "Value" or "cue" or "type".
The "attrib" is the content under the tag: eg. the string for the tag "filename"

# Run the executable
In []: object.call()

# Look through the parameter tree in python
In []: object.get_ele(keys=['key1','key2',...])

This will return the current value of the parameter in the XML associated with the path
"key1/key2/.../keyfinal/". If you use an additional input "opts='All'", then get_ele will
return all parameters along that path, giving a broader view of the XML tree.

In []: object.get_obs()
Returns a dict containing the arrays for DM, RM, and Sync (I Q U PI). If multiple sync frequencies
are used, separate arrays will be contained in a nested dict:

Ex:
array = object.get_obs()
array['Sync']['1.4']['I'] --> this will return the I array for the 1.4 GHz run. If only one freq is given,
array['Sync']['Q'] -->  then (I Q U PI) will be contained directly in the 'Sync' dict and this syntax will work.

In []: object.cleanup()
Removes entire working directory as specified in 'working_directory' at init.

****** warning: no checks are made to avoid removing the cwd or other files not related to hammurabi
'''

import os
import subprocess
import healpy as hp
import xml.etree.ElementTree as et
import random as rnd
import numpy as np
import tempfile as tf

from shutil import copyfile

class Hampy(object):


	def __init__(self,custom_parameters=None,exe_type='hamx', 
			working_directory='./test.out',hamx_path=None):

		# we offer a freedom for different executables
		if type(exe_type) is str:
			self.target = exe_type
		else:
			self.target = 'hamx'
		# find the executable
		if hamx_path is None:
			self.src_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
			executable_path = os.path.abspath('../build')
			self.executable = os.path.join(executable_path,self.target)
		else:
			self.executable = os.path.abspath(hamx_path)
		
		# create working directory and fill with a copy of the input parameter file
		self.dir = os.path.abspath(working_directory)
		if not os.path.isdir(self.dir):
			os.mkdir(self.dir)

		# nail down parameter file
		# if no input / non-string input is given, assume 'defaults.xml'
		if custom_parameters is None:
			self.base_file = os.path.join('.','defaults.xml')
		elif type(custom_parameters) is str:
			self.base_file = os.path.join('.',custom_parameters)
		else:
			print 'Error: Invalid input. Setting parameters to default file (defaults.xml)'
			self.base_file = os.path.join('.','defaults.xml')

		# copy parameter file and put this copy into working directory
		copyfile(self.base_file,working_directory + '/params.xml')
		self.temp_file = os.path.join(working_directory,'params.xml')

		# simulation output
		self.sim_map_name = {}
		self.sim_map = {}
		
		#read the copy of the parameter file into element tree
		tree = et.parse(self.temp_file)
		self.tree = tree
		self.root = self.tree.getroot()		

		self.last_call_log = ""
		self.last_call_err = ""


	def call(self,keychain=None,custom_parameter_file=None,param_path=None):
		
		import sys
		import time
		
		if keychain is not None:
			for keys in keychain:
				self.mod_par(keylist=keys)

		# replace parameter file at runtime
		self.tree.write(self.temp_file)
		if type(custom_parameter_file) is str:
			if param_path is None:	# param path is only used if specified: otherwise, assumes that param file is in cwd
				self.call_file = os.path.join('.',custom_parameter_file)
			else: 
				self.call_file = os.path.join(param_path,custom_parameter_file)
			runfile = self.call_file
		else:
			runfile = self.temp_file	
		# run hamx
		temp_process = subprocess.Popen([self.executable,runfile],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

		while temp_process.poll() is None:
			sys.stdout.write('.')
			sys.stdout.flush()
			time.sleep(0.5)
		print '\n'
		temp_process.wait()

		#if temp_process.returncode != 0:
		#	last_call_log,last_call_err = temp_process.communicate()
		#	print last_call_log
		#	print last_call_err
		
		self.last_call_log,self.last_call_err = temp_process.communicate()
		self.last_call_log=self.last_call_log.decode("utf-8")
        	if self.last_call_err: 
			self.last_call_err=self.last_call_err.decode("utf-8")
        	self.print_log(logfile='ham.log')
        	if temp_process.returncode != 0: 
			self.print_log(errfile='err.log')
        	return temp_process.returncode

		tree = et.parse(self.temp_file)
		self.tree = tree
		self.root = self.tree.getroot()
        
        def get_ele(self, keys=None, opts=None):

		root = self.root
            
                # print top parameter level if no input is given
                if keys is None or not keys or keys == ['']:
                	for child in root:
                        	print child.tag, child.attrib
                			
                else:
                   
		    if opts is None: # print the selected level
                   
                        path_str = '.' 
    		        for key in keys:
		        	path_str += '/' + key
    		
		        target = root.find(path_str)
                        print target.tag, target.attrib
                        for child in target:
                            print '|-->', child.tag, child.attrib

                    elif opts is 'All': # print all levels down to the selected level

                        path_str = '.'
			nn = 0 
    		        for key in keys:
			    nn += 1
			    path_str += '/' + key
                            target = root.find(path_str)
                            print target.tag
                            for child in target:
                                print '|'+nn*'-'+'->', child.tag, child.attrib                        
                        
                    else:
                        print ('Error: option is not supported')


	def mod_par(self, keylist=None,slct=None):
	
		if len(keylist) is 2:
			keys = keylist[0]
			tag = None
			attrib = keylist[1]
		elif len(keylist) is 3:
			keys = keylist[0]
			tag = keylist[1]
			attrib = keylist[2]
		else:
			print ('Error: invalid input')
			exit(1)
                
                if (keys is None or attrib is None):
                	print ('Error: missing input')
                	exit(1)
                
                if tag is None:
                	tag = 'value'
		
		root = self.root
		path_str = '.'
		for key in keys:
			path_str += '/' + key
    		
		target = root.find(path_str)
		if target is None:
			print ('Error: wrong element path:')
			print path_str
			exit(1)
	
		if slct is not None:
			if len(slct) is 2:
				target = root.find((path_str + "[@" + str(slct[0]) + "='" + str(slct[1]) + "']"))
			else:
				print ("Error: Invalid input -- slct must have the form ['attrib','value']")

		target.set(tag,attrib)

	
	def get_obs(self):
        # locate the files

		root = self.root
		self._get_sim_map()
		iqu_paths = {}
		self.sim_map['Sync'] = {}

		for paths in self.sim_map_name['Sync']:
			iqu_paths[paths] = os.path.join('.',self.sim_map_name['Sync'].get(paths))
		fd_path = os.path.join('.',self.sim_map_name['fd'])
		dm_path = os.path.join('.',self.sim_map_name['dm'])
		
		if(os.path.isfile(dm_path)):
			[DM] = self._read_fits_file(dm_path)
			self.sim_map['DM'] = DM

		if len(self.sim_map_name['Sync']) > 1:
			for items in iqu_paths:
				if(os.path.isfile(iqu_paths.get(items))):
					[Is,Qs,Us] = self._read_fits_file(iqu_paths.get(items))
					self.sim_map['Sync'][items] = {}
            				self.sim_map['Sync'][items]['I'] = Is
					self.sim_map['Sync'][items]['Q'] = Qs
					self.sim_map['Sync'][items]['U'] = Us
					self.sim_map['Sync'][items]['P'] = np.sqrt(np.square(Qs) + np.square(Us))
					self.sim_map['Sync'][items]['PA'] = np.arctan2(-Us,Qs)/2
		else:
			tmp_name = iqu_paths.items()
			if(os.path.isfile(tmp_name[0][1])):
				[Is,Qs,Us] = self._read_fits_file(tmp_name[0][1])
				self.sim_map['Sync']['I'] = Is
				self.sim_map['Sync']['Q'] = Qs
				self.sim_map['Sync']['U'] = Us
				self.sim_map['Sync']['P'] = np.sqrt(np.square(Qs) + np.square(Us))
				self.sim_map['Sync']['PA'] = np.arctan2(-Us,Qs)/2


		if(os.path.isfile(fd_path)):
			[Fd] = self._read_fits_file(fd_path)
			self.sim_map['Fd'] = Fd

		return self.sim_map	
	
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
	
			
	def print_log(self,start=0,stop=-1,logfile=None,errfile=None):
		'''
        	Dump the log into a file.
        	'''
        	if logfile is not None:
	        	f=open(os.path.join(self.dir,logfile),'w')
       			f.write(self.last_call_log)
        		return
        	lines=self.last_call_log.split('\n')
        	if stop==-1:
			stop=len(lines)
        	for l in range(start,stop):
			print(lines[l])

        	if errfile is not None and self.last_call_err is not None:
          		f=open(os.path.join(self.dir,errfile),'w')
			f.write(self.last_call_err)
	
	
	def cleanup(self):

		from shutil import rmtree
		# remove entire working directory	
		if os.path.isdir(self.dir):
			rmtree(self.dir)

	
	def _get_sim_map(self):
	
		root = self.root

		numsync =  len(root.findall("./Obsout/Sync[@cue='1']"))
		#print root.findall("./Obsout/Sync[@cue='1']")
#		if numsync == 1:
#			iqu_name = root.find("./Obsout/Sync[@cue='1']")
#			self.sim_map_name['iqu'] = iqu_name.get('filename')
#		else:
		self.sim_map_name['Sync'] = {}
		for sync in root.findall("./Obsout/Sync[@cue='1']"):
			self.sim_map_name['Sync'][str(sync.get('freq'))] = sync.get('filename')

		fd_name = root.find("./Obsout/Faraday")
		dm_name = root.find("./Obsout/DM")
	
		self.sim_map_name['fd'] = fd_name.get('filename')
		self.sim_map_name['dm'] = dm_name.get('filename')
		
