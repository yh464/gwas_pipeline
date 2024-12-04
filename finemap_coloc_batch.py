# -*- coding: utf-8 -*-
#!/usr/bin/env python3

'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2024-11-25

A script to batch run hyprcoloc for groups of phenotypes

Preceding workflow:
    gwa_batch.py
Required input:
    GWAS summary stats for a group of traits
    clump output (sentinel variants)
'''

def main(args):
    import os
    import pandas as pd
    
    