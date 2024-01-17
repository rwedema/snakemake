# Slurm

As a cluster workload manager, Slurm has three key functions. First, it allocates exclusive and/or non-exclusive access to resources (compute nodes) to users for some duration of time so they can perform work. Second, it provides a framework for starting, executing, and monitoring work (normally a parallel job) on the set of allocated nodes. Finally, it arbitrates contention for resources by managing a queue of pending work.&#x20;

{% embed url="https://slurm.schedmd.com/overview.html" %}

### Main Slurm Commands

* **sbatch** - submit a job script.
* **srun** - run a command on allocated compute node(s).
* **scancel** - delete a job.
* **squeue** - show state of jobs.
* **sinfo** - show state of nodes and partitions (queues).
* **smap** - show jobs, partitions and nodes in a graphical network topology.
* **scontrol** - modify jobs or show information about various aspects of the cluster

#### **sinfo**&#x20;

To find out which partitions and nodes are available and to show their status one can use `sinfo`&#x20;

```
fennaf@nuc401:~$ sinfo
PARTITION     AVAIL  TIMELIMIT  NODES  STATE NODELIST
workstations*    up 1-00:00:00      3  down* bin[206,301],nuc330
workstations*    up 1-00:00:00     99   idle bin[100-109,200-205,207-219,232-233,300,302-307,341-344],nuc[000-003,020,100-109,111-117,130,221-224,230,400-427,430]
assemblix        up 3-00:00:00      2   idle assemblix[2012,2019]
```

we can see that workstations is the default (it has a \*) and has 99 nodes idle, assemblix has only 2. Nodes are the number of computers.&#x20;

#### smap

The `smap` command is similar to the `sinfo` command, except it displays all of the information in a pseudo-graphical way.

#### squeue

The `squeue` command will report the state of running (`R`) and pending jobs (`PD`). You can use this command to find out which node your job is running on.&#x20;

```
fennaf@nuc401:~$ squeue
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
119557 assemblix snakejob   fennaf PD       0:00      4 (PartitionNodeLimit)
```

In the example above the status is `PD`, because 4 nodes are requested, but only 2 are available. In other words, it is waiting for more nodes. This will never happen, so we need to cancel this job

#### scancel

The `scancel` command will terminate pending and running job steps. It needs a `jobid`. In our example above the `jobid` is 1195557

```
fennaf@nuc401:~$ scancel 119557
fennaf@nuc401:~$ squeue
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
```

#### scontrol

We can view more details by the `scontrol` command, which is used either to tweak or display information about jobs, partition structures, and nodes.&#x20;

#### sbatch

The `sbatch` command submits a batch processing job to the slurm queue manager. These scripts typically contain one or more `srun` commands to queue jobs for processing. If you want to combine it with snakemake, you need a config file with the batch parameters like partition, nodes, tasks etc.&#x20;

```bash
snakemake --cluster "sbatch  -t {cluster.time} -N {cluster.nodes}" --cluster-config cluster_config.yaml --jobs 2 --snakefile NCBI.smk 
```

```
__default__:
  time: 01:00:00
  nodes: 4
  ntask-per-node: 10 #Request n cores be allocated per node.
  output: fenna_tries_slurm-%j.out
  error: fenna_tries_slurm-%j.error

```

#### srun

The `srun` command is used to submit jobs for execution or to initiate steps of jobs in real-time.&#x20;



### Picking the best resource for your job

Choosing the correct configuration when you submit your job will maximize its performance.&#x20;

The best way to choose a configuration is to do a test run of your job on each partition with a **small** subset of your data and benchmark it.&#x20;



### scontrol commands

#### Job Info

The `scontrol` command can be used to display information about submitted jobs, running jobs, and very recently completed jobs.&#x20;

> `scontrol show job <jobid here>`

```bash
fennaf@bin201:~$ scontrol show job 119596
JobId=119596 JobName=snakejob.download_and_count.1.sh
   UserId=fennaf(1600) GroupId=staff(1106) MCS_label=N/A
   Priority=4294900906 Nice=0 Account=(null) QOS=normal
   JobState=COMPLETED Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:42 TimeLimit=01:00:00 TimeMin=N/A
   SubmitTime=2022-03-22T14:26:52 EligibleTime=2022-03-22T14:26:52
   AccrueTime=2022-03-22T14:26:52
   StartTime=2022-03-22T14:26:52 EndTime=2022-03-22T14:27:34 Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2022-03-22T14:26:52
   Partition=workstations AllocNode:Sid=bin201:2370157
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=bin[100-103]
   BatchHost=bin100
   NumNodes=4 NumCPUs=4 NumTasks=0 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=4,node=4,billing=4
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=0 MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/homes/fennaf/Education/bfvh4dsp1/WC05/.snakemake/tmp.aqc4maq6/snakejob.download_and_count.1.sh
   WorkDir=/homes/fennaf/Education/bfvh4dsp1/WC05
   StdErr=/homes/fennaf/Education/bfvh4dsp1/WC05/slurm-119596.out
   StdIn=/dev/null
   StdOut=/homes/fennaf/Education/bfvh4dsp1/WC05/slurm-119596.out
   Power=
   NtasksPerTRES:0

```

#### Partition Info

If you would like to know detailed information about a partition, you can also use `scontrol` for that

> scontrol show partition \<partition name here>

```bash
fennaf@bin201:~$ scontrol show partition assemblix
PartitionName=assemblix
   AllowGroups=ALL AllowAccounts=ALL AllowQos=ALL
   AllocNodes=ALL Default=NO QoS=N/A
   DefaultTime=NONE DisableRootJobs=NO ExclusiveUser=NO GraceTime=0 Hidden=NO
   MaxNodes=UNLIMITED MaxTime=3-00:00:00 MinNodes=0 LLN=NO MaxCPUsPerNode=UNLIMITED
   Nodes=assemblix201[2,9]
   PriorityJobFactor=1 PriorityTier=1 RootOnly=NO ReqResv=NO OverSubscribe=NO
   OverTimeLimit=NONE PreemptMode=OFF
   State=UP TotalCPUs=240 TotalNodes=2 SelectTypeParameters=NONE
   JobDefaults=(null)
   DefMemPerNode=UNLIMITED MaxMemPerNode=UNLIMITED
```

#### Node Info

If you would like to know more about a specific node, say one of the nodes from the assemblix partition, you could get the list of nodes using the above command for getting partition info. With the list of nodes from above, put in the name of one, or more, of the nodes below to get detailed information. You can get the number of CPUs in a node, the amount of memory a node has, and the features it supports using this command.

> `scontrol show node <node name here>`

```
fennaf@bin201:~$ scontrol show node assemblix2019
NodeName=assemblix2019 Arch=x86_64 CoresPerSocket=20 
   CPUAlloc=0 CPUTot=160 CPULoad=0.06
   AvailableFeatures=assemblix
   ActiveFeatures=assemblix
   Gres=(null)
   NodeAddr=assemblix2019 NodeHostName=assemblix2019 Version=20.11.4
   OS=Linux 5.10.0-11-amd64 #1 SMP Debian 5.10.92-1 (2022-01-18) 
   RealMemory=902235 AllocMem=0 FreeMem=798473 Sockets=4 Boards=1
   CoreSpecCount=1 CPUSpecList=158-159 MemSpecLimit=1000
   State=IDLE ThreadsPerCore=2 TmpDisk=0 Weight=1 Owner=N/A MCS_label=N/A
   Partitions=assemblix 
   BootTime=2022-02-24T14:24:55 SlurmdStartTime=2022-03-18T15:26:22
   CfgTRES=cpu=160,mem=902235M,billing=160
   AllocTRES=
   CapWatts=n/a
   CurrentWatts=0 AveWatts=0
   ExtSensorsJoules=n/s ExtSensorsWatts=0 ExtSensorsTemp=n/s
   Comment=(null)

```
