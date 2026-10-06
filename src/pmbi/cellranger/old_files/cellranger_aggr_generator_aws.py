import re
import scriptgen
import os
import cellranger_command
import mbiaws

grp = cellranger_command.CellrangerAggrGroups()
for i in mbiaws.s3.lib.list_object_keys("microbioinfo-storage", "betts_t1d/032723__HWTLWDSX5__/cellranger_multi/"):
    if mbiaws.s3.lib.object_key_matches("[.]outs[.]tar[.]*(gz)*$", i):
        s = re.search(cellranger_command.super_sample_p, os.path.basename(i))
        super_sample = s.group(0)
        pid = cellranger_command.get_patient_id(super_sample)
        grp.insert(pid, i)
setup_cmds = ["source /shared-ebs/source_me.bash",
              "REPO=/shared-ebs/microbioinfo-aws/scripts/",
              "unset AWS_PROFILE; unset AWS_SHARED_CREDENTIALS_FILE", # Make sure boto3 used IAM profile creds
              "mkdir -p /local-ebs/cellranger/multi",
              "mkdir -p /local-ebs/cellranger/aggr",
              "cd /local-ebs/cellranger/multi"]
for i,g in enumerate(grp.groups.items()):
    g,k = g
    aggr = cellranger_command.CellrangerAggrCommand.from_s3_archive_keys(g,k,mode="multi")
    sg = scriptgen.SlurmScriptGenerator(
            jobname=f"craggr_{i}",
            nodes=1,
            tasks_per_node=1,
            cpus_per_task=16,
            mem=120,
            time=168,
            partition="r5ad4x"
            )

    for sc in setup_cmds:
        sg.add_command(sc)
    for ad in aggr.awsdirs():
        sg.add_command(f"python $REPO/s3-download-multi.py --bucket microbioinfo-storage --prefix {os.path.join(ad, f'{aggr.pid}_')} --chunksize 1.024e8 --pattern '[.]multi[.]outs[.]tar[.]gz$'")
    sg.add_command(f"python $REPO/jobgen/cellranger_untar_multi_outs.py {' '.join(aggr.to_fetch_bases)}")
    sg.add_command("cd /local-ebs/cellranger/aggr")
    sg.add_command(f"cat << EOF > aggr.csv\n{aggr.csv}\nEOF")
    sg.add_command(f"cellranger aggr --csv aggr.csv --id {aggr.pid} --localcores 16 --localmem 120")
    sg.add_command(f"tar cvfz {aggr.pid}.aggr.outs.tar.gz {os.path.join(aggr.pid, 'outs')}")
    sg.add_command(f"python $REPO/s3-multipart-upload.py --bucket microbioinfo-storage --key {os.path.join(aggr.awsdirs().pop().replace('cellranger_multi', 'cellranger_aggr'), f'{aggr.pid}.aggr.outs.tar.gz')} --partsize 1000000000 --nproc 16 --largefile {aggr.pid}.aggr.outs.tar.gz")
    sg.write()

    

