import os
from xmlrpclib import ServerProxy
import pathlib2 as pl
import enviroment
import analysis

def connect(ip, port):
    pymol = ServerProxy(uri = 'http://%s:%d/RPC2' % (ip, port))
    return pymol

class Pymol(object):
    def __init__(self, analytics, protein_path):
        self._ana = analytics
        self.prot_path = protein_path
        self.py_con = connect(ip='localhost')

    def setup(self):
        self.py_con.delete('all')
        self.py_con.load(self.prot_path)
        self.py_con.show_as('cartoon')
        self.py_con.color('green')

    # The counter list argument is used if you want to use a
    # pre-prepared list of contact counts like if you use the pipeline,
    # otherwise it will take the average conact counts of all runs for
    # the specific job.
    def show_bindings(self, binding = "multi", counter_list = None):
        if counter_list is not None:
            self._pymol_bindings(counter_list)
        elif binding == "single":
            self._pymol_bindings(self._ana._contact_frequencies_protein(
                self._ana.dataframe))
        elif binding == "multi":
            self._pymol_bindings(self._ana._contact_frequencies_protein(
                self._ana._multi_concat))
        else:
           raise ValueError('Only "multi" and "single" binding probability is supported')

def _pymol_bindings(self, df):
    if df.has_key(0):
        del df[0]
    gradient = None
    for residue in df.keys():
        if int(df[residue]*100) in range(1, 5):
            gradient = 'yellow'
        elif int(df[residue]*100) in range(5, 10):
            gradient = 'orange'
        elif int(df[residue]*100) in range(10, 20):
            gradient = 'brown'
        elif int(df[residue]*100) >= 20:
            gradient = 'red'
        else:
            pass
        self.py_con.select('resi %i' % residue)
        self.py_con.color(gradient, 'sele')
        self.py_con.show('spheres', 'sele')
     self.py_con.deselect()
