#!/bin/bash
#SBATCH --job-name=rclone
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=35gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=24:00:00
#SBATCH --output=rclone.out
#SBATCH --error=rclone.err

# Load required modules
module load rclone

tar cf - /orange/zhangw/NS3842 | pigz -p 16 > NS3842.tar.gz
rclone copy NS3842.tar.gz dbox:ProjectTrackingApp

