#!/bin/bash
#
# Most tools are already installed in /ebs1. Thias will install some
# missing tools with apt.
#
# Tested on the "Ubuntu Server 18.04 LTS (HVM), SSD Volume Type - ami-005bdb005fb00e791 (64-bit x86)" image
#

#Install tools with apt
sudo apt update
sudo apt upgrade -qyy
sudo apt install -y make pkg-config emacs virtualenv
sudo apt install -y samtools
sudo apt install -y tabix
sudo apt install -y libssl-dev
#sudo apt install -y docker.io
#sudo usermod -aG docker $USER
sudo apt install -y last-align
sudo apt install -y snakemake
