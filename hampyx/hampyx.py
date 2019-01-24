# -*- coding: utf-8 -*-
'''
Developed by Joe Taylor, based on the initial work by Theo Steininger and Jiaxin Wang
Reviewed by Jiaxin Wang

temporary/testing version of hammurabiX python wrapper

warning:

working directory is set as the same directory as this file
hammurabiX executable path is by default /usr/local/hammurabi/bin/hamx
it relies on subprocess to fork c++ routine

methods:

# Import class
In []: import hampyx as hpx

# Initialize object
# exe_name: the name of executable, 'hampyx' by default
In []: object = hampyx.Hampy(exe_name='')

# Modify parameter value from base xml file to temp xml file
In []: object.mod_par(keychain=['key1','key2',...],tag='',attrib='')
The strings 'key1', 'key2', etc represent the path to the desired
parameter, going through the xml.

The "tag" is the label for the parameter: eg. "Value" or "cue" or "type".
The "attrib" is the content under the tag: eg. the string for the tag "filename"

# Run the executable
In []: object.call()

# Look through the parameter tree in python
In []: object.print_par(keychain=['key1','key2',...])

This will return the current value of the parameter in the XML associated with the path
"key1/key2/.../keyfinal/". If you use an additional input "opt='all'", then get_ele will
return all parameters along that path, giving a broader view of the XML tree.

In []: object.get_obs()
Returns a dict containing the arrays for DM, RM, and Sync (I,Q,U,PI,PA).
If multiple sync frequencies are used, separate arrays will be contained in a nested dict:

Ex:
array = object.get_obs()
array['sync']['1.4']['I'] --> this will return the I array for the 1.4 GHz run.

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

class hampyx(object):

        # default executable path is '/usr/local/hammurabi/bin/hamx'
        # default executable path is './params.xml'
	def __init__(self,
                     exe_path=None,
                     xml_path=None):
                # current working directory
                self.wk_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
		# find the executable
		# by default hammurabiX executable "hamx" is in the same directory as this file
		if exe_path is None:
			self.executable = os.path.abspath('/usr/local/hammurabi/bin/hamx')
		else:
			self.executable = os.path.abspath(exe_path)#scopend
		if xml_path is None:
                        self.base_file = os.path.join(self.wk_dir,'params.xml')
                else:
                        self.base_file = os.path.abspath(xml_path)#scopend
                # assign tmp file name by base file name
                self.temp_file = self.base_file
		# read from base parameter file
		self.tree = et.parse(self.base_file)
		# simulation output
		self.sim_map_name = {}
		self.sim_map = {}	
                #
		self.last_call_log = ""
		self.last_call_err = ""
		
	def call(self,verbose=False):
		import sys
		import time
		# create new temp parameter file
		if self.temp_file is self.base_file:
                        self._new_xml_copy()
		
		temp_process = subprocess.Popen([self.executable,self.temp_file],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
		temp_process.wait()
                # if need verbose output
		if verbose is True:
                        self.last_call_log,self.last_call_err = temp_process.communicate()
                        self.last_call_log = self.last_call_log.decode("utf-8")
                        self.last_call_err = self.last_call_err.decode("utf-8")
                        self.write_log(filename='hampy.log')
                # if quiet, only print upon error
                else:
                        if temp_process.returncode != 0:
                                last_call_log,last_call_err = temp_process.communicate()
                                print last_call_log
                                print last_call_err
		# grab output maps and delete temp files
		self._get_sims()
		self._del_xml_copy()
		
	# grab simulation output from disk and delete corresponding fits file
	def _get_sims(self):
        	# locate the files
        	dm_path = self.sim_map_name['dm']
        	fd_path = self.sim_map_name['fd']
        	sync_paths = {}
                # sync_path has multiple contents
                for i in self.sim_map_name['sync']:
                        sync_paths[i] = self.sim_map_name['sync'].get(i)
                # read dispersion measure and delete file
		if(os.path.isfile(dm_path)):
			[DM] = self._read_fits_file(dm_path)
			self.sim_map['DM'] = DM
			os.remove(dm_path)
		else:
                        print 'ERR: missing file'
                        exit(1)
                # read faraday depth and delete file
		if(os.path.isfile(fd_path)):
			[Fd] = self._read_fits_file(fd_path)
			self.sim_map['F'] = Fd
          		os.remove(fd_path)
          	else:
                        print 'ERR: missing file'
                        exit(1)
		# read synchrotron pol. and delete file
                for i in sync_paths:
                        # if file exists
                        if(os.path.isfile(sync_paths.get(i))):
                                [Is,Qs,Us] = self._read_fits_file(sync_paths(i))
                                # build top level map for nesting
                                self.sync_map['sync'][i] = {}
                                self.sync_map['sync'][i]['I'] = Is
                                self.sync_map['sync'][i]['Q'] = Qs
                                self.sync_map['sync'][i]['U'] = Us
                                # polarisation intensity
                                self.sync_map['Sync'][items]['PI'] = np.sqrt(np.square(Qs) + np.square(Us))
                                # polarisatioin angle, IAU convention
                                self.sync_map['Sync'][items]['PA'] = np.arctan2(Us,Qs)/2.0
                                os.remove(sync_paths(i))
                        else:
                                print 'ERR: missing file'
                                exit(1)
		if(os.path.isfile(iqu_path)):
			[Is,Qs,Us] = self._read_fits_file(iqu_path)
            		self.sim_map['I'] = Is
			self.sim_map['Q'] = Qs
			self.sim_map['U'] = Us
			os.remove(iqu_path)
	
	# read a single fits file with healpy
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

        # make a temporary parameter file copy and rename output file with random mark
	def _new_xml_copy(self):
		# create a random file name which doesn't exist currently
		fd,new_path = tf.mkstemp(prefix='params_',suffix='.xml',dir=self.wk_dir)
		os.close(fd)
		rnd_idx = new_path[new_path.index('params_')+7:-4]
		self.temp_file = new_path
		# copy base_file to temp_file
		root = self.tree.getroot()
		# count number of sync output files
		syncnum = len(root.findall("./Obsout/Sync[@cue='1']"))
		# for each sync output
		for sync in root.findall("./Obsout/Sync[@cue='1']"):
                        freq = str(sync.get('freq'))
                        self.sim_map_name['sync'][freq] = os.path.join(self.wk_dir,'iqu_sync_'+freq+'_'+rnd_idx+'.fits')
                        sync.set('filename',self.sim_map_name['sync'][freq])#scopend
		self.sim_map_name['fd'] = os.path.join(self.wk_dir,'fd_'+rnd_idx+'.fits')
		self.sim_map_name['dm'] = os.path.join(self.wk_dir,'dm_'+rnd_idx+'.fits')
		root.find('./Obsout/Faraday').set('filename',self.sim_map_name['fd'])
		root.find('./Obsout/DM').set('filename',self.sim_map_name['dm'])
		# automatically create a new file
		tree.write(self.temp_file)
	
        # delete temporary parameter file copy
	def _del_xml_copy(self):
		if self.temp_file is self.base_file:
			print 'ERR: _del_xml_copy failed'
			exit(1)
		else:
			os.remove(self.temp_file)
			self.temp_file = self.base_file

        # print log file, append to existed file instead of rewriting
        def write_log(self,filename=None):
                f = open(os.path.join(self.wk_dir,filename),'a')
                f.write('CALL_LOG:\n')   
                f.write(self.last_call_log)
                f.write('ERR_LOG:\n')
                f.write(self.last_call_err)

        # modify parameter in self.tree
	def mod_par(self,tanent=None):
                if len(tanent) is 2:
                        keys = tanent[0]
                        tag = 'value'
                        attrib = tanent[1]#scopend
                elif len(tanent) is 3:
                        keys = tanent[0]
                        tag = tanent[1]
                        attrib = tanent[2]#scopend
                else:
			print 'ERR: missing input'
			exit(1)#scopend
		root = self.tree.getroot()
		path_str = '.'
		for key in keys:
			path_str += '/' + key#scoped
		target = root.find(path_str)
		if target is None:
			print 'ERR: wrong element path:'
			print path_str
			exit(1)#scopend
		target.set(tag,attrib)
		
        # print a certain parameter
        # argument of type ['path','to','key']
	def print_par(self,keychain=None,opt=None):
                root = self.tree.getroot()
                # print top parameter level if no input is given
                if keychain is None:
                        for child in root:
                                print child.tag, child.attrib
                else:
                    if opt is None: # print the selected level
                        path_str = '.' 
    		        for key in keychain:
		        	path_str += '/' + key#scopend
		        target = root.find(path_str)
                        print target.tag, target.attrib
                        for child in target:
                            print '|-->', child.tag, child.attrib
                    elif opt is 'all': # print all levels down to the selected level
                        path_str = '.' 
    		        for key in keychain:
			    path_str += '/' + key
                            target = root.find(path_str)
                            print target.tag
                            for child in target:
                                print ('|-->', child.tag, child.attrib)
                    else:
                        print ('ERR: unsupported option')
        
        # deletes an parameter and all of its children
        # USE WITH CAUTION
	# argument is the path to the element (e.g. ['Grid','SunPosition','x'])
	def del_par(self, keychain=None):
		root = self.tree.getroot()
		if keys is not None:
			path_str = '.'
			n = 1
			for key in keychain:
				if n is len(keychain):
                                        par_path_str = path_str#scopend
				path_str += '/' + key
				n += 1#scopend
			target = root.find(path_str)
			parent = root.find(par_path_str)
			parent.remove(target)
		else:
			print ('ERR: missing input')
