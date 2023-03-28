#! /Users/mshirota/opt/anaconda3/bin/python3

import sys, os, re, types, gzip
from argparse import ArgumentParser
#import pandas as pd
#import numpy as np

ss3dict = { 'H': 'H', 'G': 'H', 'I': 'H', 'B': 'S', 'E' : 'S', ' ' : 'L', 'S' : 'L', 'T' : 'L' }
maxasadict = { 'A' : 115, 'R' : 225, 'N' : 160, 'D' : 150, 'C' : 135, 'Q' : 180, 'E' : 190, 'G' : 75, 'H': 195, 'I': 175,\
               'L' : 170, 'K' : 200, 'M' : 185, 'F' : 210, 'P' : 145, 'S' : 115, 'T' : 140, 'W' : 255, 'Y': 230, 'V':155 }

header = [ 'Chain', 'ResID', 'AA', 'SS', 'SS3', 'ASA', 'rASA', 'Phi','Psi', 'ABEG' ]

def read_dssp( f ):

    start = False
    data = []
    for l in f:
        if l.startswith('  #  RESIDUE AA STRUCTURE BP1 BP2  ACC'):
            start = True
            continue
        if not start:
            continue
        resid, chain, aa, ss, acc, phi, psi = l[5:10], l[11:12], l[13:14], l[16:17], l[35:38], l[103:109], l[109:115]
        if aa == '!':
            continue
        if re.match( '^[a-z]$', aa ):
            aa = 'C'
        ss3 = ss3dict.get( ss, 'L' )
        rasa = '{:.3f}'.format( float( acc ) / maxasadict[ aa ] )
        abeg = 'A' if float(psi) >= -75 and float(psi) < 50 else 'B' if float(phi) < 0 \
            else 'G' if float(psi) >= -100 and float(psi) < 100 else 'E'
        data.append( { k:v for k, v in zip( header, [ chain, resid, aa, ss, ss3, acc, rasa, phi, psi, abeg ] ) } )
    return data

def make_dssp_dict( data ):
    dssp_dict = {}
    for d in data:
        dssp_dict[ ( d['Chain'], d['ResID'] ) ] = d
    return dssp_dict

def calc_dasa( h, d, md ):

    if md is None:
        return None
    if h == 'drASA':
        return '{:.3f}'.format( float( md[h[1:]] ) - float( d[h[1:]] ) )
    else:
        return '{}'.format( int( md[h[1:]] ) - int( d[h[1:]] ) )

def _main():
    parser = ArgumentParser()
    parser.add_argument( '-i', '--ifilename', dest='ifilename', help='input filename, default=', default='')
    parser.add_argument( '-m', '--monomerfilename', dest='monomerfilename', help='input filename, default=None', default=None)
    parser.add_argument( '-o', '--ofilename', dest='ofilename', help='output filename, default=', default='')
    #parser.add_argument( '-I', '--idirname', dest='idirname', help='input directory name, default=./', default='')
    #parser.add_argument( '-O', '--odirname', dest='odirname', help='output directory name, default=output/', default='output/')
    parser.add_argument( '-c', '--chain', dest='chain', help='chain ids that are  shown in the output, default=None (all chains)', default=None )
    parser.add_argument( '-r', '--residue', dest='residue', help='residue ids that are shown in the output, default=None (all residue ids)', default=None )
    parser.add_argument( '-f', '--feature', dest='feature', help='feature columns that are shown in the output, default=None (all columns)', default=None )


    args = parser.parse_args()

    for attr in dir( args ):
        if attr.startswith('_'):
            continue
        if type( getattr( args, attr ) ) == str and  not re.match('^[a-zA-Z0-9._/\,\-]*$', getattr( args, attr ) ):
            raise ValueError('Bad Input for Argument:{attr}'.format( attr = attr ))

    if args.monomerfilename is not None:
        with gzip.open( args.monomerfilename, 'rt' ) if args.monomerfilename.endswith('.gz') else open( args.monomerfilename, 'rt' ) as mf:
            mdata = make_dssp_dict( read_dssp( mf ) )
        outheader = header + [ 'dASA', 'drASA']
    else:
        if args.feature is not None and ( 'dASA' in args.feature.split(',') or 'drASA' in args.feature.split(',') ):
            print( 'Error, cannot calculate dASA or drASA unless monomerfile is given!', file = sys.stderr )
            sys.exit(1)
        mdata = None
        outheader = header
        

        
    if args.feature is not None:
        for feature in args.feature.split(','):
            if feature not in outheader:
                print( f'feature {feature} is not allowed!', file = sys.stderr )
                sys.exit(1)
        
    with gzip.open( args.ifilename, 'rt' ) if args.ifilename.endswith('.gz') else open( args.ifilename, 'rt' ) if args.ifilename and args.ifilename != '-' else sys.stdin as f, gzip.open( args.ofilename, 'wt' ) if args.ofilename.endswith('.gz') else open( args.ofilename, 'wt' ) if args.ofilename and args.ofilename != '-' else sys.stdout as of:
        print( '\t'.join( [ h for h in outheader if args.feature is None or h in ( [ 'Chain', 'ResID', 'AA' ] + args.feature.split(',') ) ] ), file = of )
        data = read_dssp( f )
        for d in data:
            if ( args.residue is None or d['ResID'].strip() in args.residue.split(',')  ) and ( args.chain is None or d['Chain'] in args.chain.split(',') ):
                print( '\t'.join( [ d[h] if h not in ['drASA', 'dASA'] else calc_dasa( h, d, mdata.get( ( d['Chain'], d['ResID'] ), None ) ) for h in outheader if args.feature is None or h in ( [ 'Chain', 'ResID', 'AA' ] + args.feature.split(',') ) ] ), file = of )

if __name__ == '__main__':
    _main()

