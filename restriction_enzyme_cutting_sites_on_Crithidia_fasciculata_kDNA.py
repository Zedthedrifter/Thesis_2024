#!/usr/bin/env python3


#%% MODULES TO IMPORT 

from __future__ import division
import matplotlib.pyplot as plt
import os
import re
import sys
import numpy
import argparse
import subprocess as sbp
import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

kdna='/home/zed/disk1/Crithidia/minicircle_assembly/Cf.kDNA.fasta'
kdna=SeqIO.to_dict(SeqIO.parse(kdna,'fasta'))
#kdna['maxi']=kdna['Cf_k99_contigk129_1263_len24398']
#kdna.pop('Cf_k99_contigk129_1263_len24398')


def make_ezm_dict():
  ezm_dict={}
  reseq=open('enzyme.txt')
  for l in reseq:
    l=l.strip('\n').split('\t')
    ezm_dict[l[1]]={}
    if l[0].startswith('('):
      ezm_dict[l[1]]['seq']=l[0].replace('/','').split(')')[1].split('(')[0]
      print(l[1],l[0].replace('/','').split(')')[1].split('(')[0])
    else:
      ezm_dict[l[1]]['seq']=l[0].replace('/','').split('(')[0]
  return(ezm_dict)

ambiguous={'N':'[GACT]',
'R':'[AG]',
'Y'	:'[CT]',
'S':	'[GC]',
'W':'[AT]',
'K':'[GT]',
'M':'[AC]',
'B':'[CGT]',
'D':'[AGT]',
'H':'[ACT]',
'V':'[ACG]'}

def find_sites(ezm_dict):
  for ezm in ezm_dict:
    #if re.search(r'[^AGCT]',ezm_dict[ezm]['seq']):
    #  print(f"ambiguious bases :{ezm}\t{ezm_dict[ezm]['seq']}")
    motif=''.join([ambiguous.get(b,b) for b in ezm_dict[ezm]['seq']])
    #print(motif)
    hits=[]
    for m in kdna:
      if re.search(motif,str(kdna[m].seq)) or re.search(motif,str(kdna[m].seq.reverse_complement())):
        hits.append(m)
    ezm_dict[ezm]['hits']=hits
  return(ezm_dict)


def find_sites_details(outfile='restriction_site_details.csv'):
  site_dict,order={},['seq']+list(kdna.keys())
  for ezm in ezm_dict:
    #if re.search(r'[^AGCT]',ezm_dict[ezm]['seq']):
      #print(f"ambiguious bases :{ezm}\t{ezm_dict[ezm]['seq']}")
    motif=''.join([ambiguous.get(b,b) for b in ezm_dict[ezm]['seq']])
    site_dict[ezm]={k:'' for k in order}
    site_dict[ezm]['seq']=ezm_dict[ezm]['seq']
    hits=[]
    for m in kdna:
      if re.search(motif,str(kdna[m].seq)) or re.search(motif,str(kdna[m].seq.reverse_complement())):
        site_dict[ezm][m]={'forward':[],'rc':[]}
        for match in re.finditer(motif,str(kdna[m].seq)):
          site_dict[ezm][m]['forward'].append(match.start())
        for match in re.finditer(motif,str(kdna[m].seq.reverse_complement())):
          site_dict[ezm][m]['rc'].append(match.start())
  df=pd.DataFrame.from_dict(site_dict,orient='index')
  df.to_csv(outfile)
  #print(order)
  return(site_dict)
 
## + and -
def find_sites_plus_minus(outfile='restriction_site_plus_minus.csv'):
  site_dict,order={},['seq']+list(kdna.keys())
  for ezm in ezm_dict:
    motif=''.join([ambiguous.get(b,b) for b in ezm_dict[ezm]['seq']])
    site_dict[ezm]={k:'' for k in order}
    site_dict[ezm]['seq']=ezm_dict[ezm]['seq']
    hits=[]
    for m in kdna:
      if re.search(motif,str(kdna[m].seq)) or re.search(motif,str(kdna[m].seq.reverse_complement())):
        site_dict[ezm][m]='+'
      else:
        site_dict[ezm][m]='-'
  df=pd.DataFrame.from_dict(site_dict,orient='index')
  df.to_csv(outfile)
  #print(order)
  return(site_dict)
#
#
#
def output_hits(ezm_dict,site_dict,output='restriction_site_hits.txt'):
  f=open(output,'w')
  #minicircles only
  f.write('Enzymes that cut minicircles only:\n\n')
  for k in ezm_dict:
    if 'maxi' not in ezm_dict[k]['hits'] and ezm_dict[k]['hits'] !=[]:
      tmp={m: {'f':len(site_dict[k][m]['forward']),'r':len(site_dict[k][m]['rc'])} for m in ezm_dict[k]['hits']}
      f.write(f"{k} \t{ezm_dict[k]['seq']}\ttargeted circles({len(ezm_dict[k]['hits'])}):\n")
      for m in tmp:
        f.write(f"\t{m}:'forward':{tmp[m]['f']}\t'reverse complement':{tmp[m]['r']}\n")
  #maxicircle only
  f.write('\nEnzymes that cut only maxicircles:\n\n')
  for k in ezm_dict:
    if ezm_dict[k]['hits']==['maxi']:
      tmp={m: {'f':len(site_dict[k][m]['forward']),'r':len(site_dict[k][m]['rc'])} for m in ezm_dict[k]['hits']}
      f.write(f"{k} \t{ezm_dict[k]['seq']}\ttargeted circles({len(ezm_dict[k]['hits'])}):\n")
      for m in tmp:
        f.write(f"\t{m}:'forward':{tmp[m]['f']}\t'reverse complement':{tmp[m]['r']}\n")
  #maxi and mini
  f.write('\nEnzymes that cut both maxi and minicircles:\n\n')
  for k in ezm_dict:
    if 'maxi' in ezm_dict[k]['hits'] and ezm_dict[k]['hits']!=['maxi']:
      tmp={m: {'f':len(site_dict[k][m]['forward']),'r':len(site_dict[k][m]['rc'])} for m in ezm_dict[k]['hits']}
      f.write(f"{k} \t{ezm_dict[k]['seq']}\ttargeted circles({len(ezm_dict[k]['hits'])}):\n")
      for m in tmp:
        f.write(f"\t{m}:'forward':{tmp[m]['f']}\t'reverse complement':{tmp[m]['r']}\n")
      #f.write(f"{k} \t{ezm_dict[k]['seq']}\ttargeted circles:{ezm_dict[k]['hits']}\n")
  #no cutting sites
  f.write('\nnot cutting site is found \n\n')
  for k in ezm_dict:
    if ezm_dict[k]['hits'] ==[]:
       f.write(f"{k}\t{ezm_dict[k]['seq']}\n")
  f.close()
 
     
ezm_dict=make_ezm_dict()
ezm_dict=find_sites(ezm_dict)
#site_dict=find_sites_details()
#output_hits(ezm_dict,site_dict)
site_dict=find_sites_plus_minus(outfile='restriction_site_plus_minus.csv')