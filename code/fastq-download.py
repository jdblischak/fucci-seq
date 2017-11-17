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
#   python code/fastq-download.py data/md5/YG-PYT-Fucci1.md5 Genomics/NGS-2017/171109_700819F_0583_ACAWRYACXX-YG-PYT-Fucci1 /project2/gilad/fucci-seq/fastq

import docopt
import getpass
import hashlib
import os
import pysftp
import sys

def main(md5file, remotedir, outdir = ".", hostname = "fgfftp.uchicago.edu", username = "gilad"):

    # Connect to server
    p = getpass.getpass("Authenticate: ")
    sftp = pysftp.Connection(hostname, username = username, password = p)

    # Import md5 checksums computed by core
    md5_core = open(md5file, "r")

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    for line in md5_core:
        cols = line.strip().split()
        md5 = cols[0]
        fname = os.path.basename(cols[1])
        if "Undetermined" in fname or fname[-8:] != "fastq.gz":
            continue
        sys.stdout.write("Downloading %s\n"%(fname))
        # Organize the FASTQ files into subdirectories based on the C1 chip
        chip = fname.split("-")[3]
        outdir = outdir + "/" + chip
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        localpath = outdir + "/" + fname
        remotepath = remotedir + "/" + fname
        if sftp.exists(remotepath):
            sftp.get(remotepath, localpath)
        else:
            sys.stderr.write("Does not exist on remote server:\t%s\n"%(remotepath))
            continue
        if not os.path.exists(localpath):
            sys.stderr.write("Download failed:\t%s\n"%(remotepath))
            continue
        with open(localpath, "rb") as fq:
            md5_local = hashlib.md5(fq.read()).hexdigest()
            if md5_local != md5:
                sys.stderr.write("Download incomplete:\t%s\n"%(remotepath))

        md5_core.close()
        sftp.close()

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    print(args)
    main(md5file = args["<md5file>"],
         remotedir = args["<remotedir>"],
         outdir = args["<outdir>"])
