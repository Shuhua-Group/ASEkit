import sys
from multiprocessing import Pool
import argparse
import os 
description = "Description:\n\n" + \
              "ASEkit allows you to run convention analysis , \n" \
              "Call ASE from WASP filter bam \n" \
              "aseQTL to find cis-candidate regulatory SNP\n" \
              "Detcte an assocaition between the value of allelic imbalance of exonic SNP and traits or environment factors\n" \
              "For further information, see the help of each subcommand."

parser = argparse.ArgumentParser(description=description, prog='ASEkit')
##version
parser.add_argument('-v','--version', action='version', version='%(prog)s 1.0.1')
version='1.0.0'
def help_info():
    help_info='ASEkit : A python package apply for analyzing ASE,\n\n' + \
            'Function1: Calling'+'\n'+'Call ASE from RNAseq data and vcf file,\n\n' + \
            'Function2: Filter '+'\n'+'Filter ASE that satisfy statstic condition,\n\n' + \
            'Function3:  aseQTL '+'\n'+'Detect cis-candidate regulatory SNP,\n\n' + \
            'Function4: association '+'\n'+'Detcte an assocaition between the value of allelic imbalance of exonic SNP and traits or environment factors'
    return help_info

def main():
    if len(sys.argv)==1:
        print(help_info())
        sys.exit(1)

    try:
        
        # command for testing
        if sys.argv[1] == 'Calling':
            from  ASEkit import aseCalling
            aseCalling.main()
        elif sys.argv[1] == 'Filter':
            from  ASEkit import aseFilter
            aseFilter.main()
        elif sys.argv[1] == 'aseQTL':
            from  ASEkit import aseQTL
            aseQTL.main()
        elif sys.argv[1] == 'aseTrait':
            from ASEkit import  aseTrait
            aseTrait.main()
        elif sys.argv[1] == 'test':
            test_filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'example.data')
            print('test data filepath is: '+test_filepath)
        elif sys.argv[1] in ['version', '--version','-V', '-v']:
            print('[ASEkit] version: %s' % version)
        # command for help information
        elif sys.argv[1] in ['help', '--help', '-h']:
            print(help_info())
        else:
            print('Can not know this parameter,please see README of the package' )
    except Exception :
        print('Can not know this parameter,please see  README of the package' )
if __name__ == '__main__':
    main()
