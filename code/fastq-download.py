#!/usr/bin/env python

"""Download FASTQ files from core FTP server.

Usage:
  fastq-download.py <md5file> <remotedir> <outdir>

Arguments:
  md5file       The file with md5 checksums provided by the core.
  remotedir     The path on the remote filesystem to the FASTQ files.
  outdir        The path on the local filesystem to save the FASTQ files.

Options:
  -h --help     Show this screen.

"""

# Example usage:
#   python code/fastq-download.py data/md5/YG-PYT-Fucci1.md5 Genomics/NGS-2017/171109_700819F_0583_ACAWRYACXX-YG-PYT-Fucci1/FastQ /project2/gilad/fucci-seq/fastq
#
# Implementation details:
#
# The md5 file has 2 columns, no header, and is space-delimted. The first column
# is the md5 checksum. The second column is the file, though the path is not the
# same as on the FTP server, and thus must be discarded. Also the file contains
# the checksums for other miscellaneous files that need to be skipped. This
# script only donwloads fastq.gz files, and also exlcudes the Undetermined
# files. Here are two example lines from an input file:
#
# aed3c3f561e6e28fd16174620f13f127  /media/data1/NewSequencerRuns/171109_700819F_0583_ACAWRYACXX/Unaligned_YG-PYT-Fucci1/YG-PYT-Fucci1-20170905-B08_S20_L001_R1_001.fastq.gz
# 214b05c57da97e00e0f4a2cf2ec47767  /media/data1/NewSequencerRuns/171109_700819F_0583_ACAWRYACXX/Unaligned_YG-PYT-Fucci1/YG-PYT-Fucci1-20170905-B09_S21_L001_R1_001.fastq.gz
#
# The downloaded files are further organized on the local filesystem. Each FASTQ
# is put in a subdirectory corresponding to the C1 chip the single cell came
# from. It will look something like this:
#
# outdir
# ├── 20170905
# │   ├── YG-PYT-Fucci1-20170905-A01_S1_L001_R1_001.fastq.gz
# │   ├── YG-PYT-Fucci1-20170905-A02_S2_L001_R1_001.fastq.gz
# │   ├── YG-PYT-Fucci1-20170905-A03_S3_L001_R1_001.fastq.gz
# │   ├── YG-PYT-Fucci1-20170905-A04_S4_L001_R1_001.fastq.gz
# │   ├── YG-PYT-Fucci1-20170905-A05_S5_L001_R1_001.fastq.gz
# │   ├── YG-PYT-Fucci1-20170905-A06_S6_L001_R1_001.fastq.gz
# │   ├── YG-PYT-Fucci1-20170905-A07_S7_L001_R1_001.fastq.gz
# │   ├── YG-PYT-Fucci1-20170905-A08_S8_L001_R1_001.fastq.gz
# │   └── YG-PYT-Fucci1-20170905-A09_S9_L001_R1_001.fastq.gz
# ├── 20170906
# ├── 20170907
# ├── 20170908
# └── 20170910
#
# The script will skip files that have already been downloaded and removes those
# whose md5 checksum does not match.

import docopt
import getpass
import hashlib
import os
import pysftp
import sys

def main(md5file, remotedir, outdir = ".", hostname = "fgfftp.uchicago.edu",
         username = "gilad"):

    # Connect to server
    p = getpass.getpass("Authenticate: ")
    sftp = pysftp.Connection(hostname, username = username, password = p)

    # Import md5 checksums computed by core
    md5_core = open(md5file, "r")

    # Download each file individually and verify the md5 checksum
    for line in md5_core:
        cols = line.strip().split()
        md5 = cols[0]
        fname = os.path.basename(cols[1])
        if "Undetermined" in fname or fname[-8:] != "fastq.gz":
            continue
        sys.stdout.write("Downloading %s\n"%(fname))
        # Organize the FASTQ files into subdirectories based on the C1 chip
        chip = fname.split("-")[3]
        outdir_chip = outdir + "/" + chip
        os.makedirs(outdir_chip, exist_ok = True)
        localpath = outdir_chip + "/" + fname
        remotepath = remotedir + "/" + fname
        #import ipdb; ipdb.set_trace()
        if not sftp.exists(remotepath):
            sys.stderr.write("Does not exist on remote server:\t%s\n"%(remotepath))
            continue
        elif os.path.exists(localpath):
            sys.stderr.write("Already exists:\t%s\n"%(remotepath))
            continue
        else:
            sftp.get(remotepath, localpath)
        if not os.path.exists(localpath):
            sys.stderr.write("Download failed:\t%s\n"%(remotepath))
            continue
        with open(localpath, "rb") as fq:
            md5_local = hashlib.md5(fq.read()).hexdigest()
            if md5_local != md5:
                sys.stderr.write("Download incomplete:\t%s\n"%(remotepath))
                os.remove(localpath)

    md5_core.close()
    sftp.close()

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    main(md5file = args["<md5file>"],
         remotedir = args["<remotedir>"],
         outdir = args["<outdir>"])
