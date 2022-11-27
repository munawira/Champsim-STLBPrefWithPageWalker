<p align="center">
  <h1 align="center"> ChampSim </h1>
  <p> ChampSim is a trace-based simulator for a microarchitecture study. You can sign up to the public mailing list by sending an empty mail to champsim+subscribe@googlegroups.com. If you have questions about how to use ChampSim, you can often receive a quicker response on the mailing list. Please reserve GitHub Issues for bugs. <p>
</p>

# To Clone Baseline ChampSim repository
```
git clone https://github.com/ChampSim/ChampSim.git
```

#To clone this project

```
git clone https://github.com/munawira/Champsim-STLBPrefWithPageWalker
```


# Compile

ChampSim takes a JSON configuration script. Examine `champsim_config.json` for a fully-specified example. All options described in this file are optional and will be replaced with defaults if not specified. The configuration script can also be run without input, in which case an empty file is assumed.
For this project, if you want to use a prefetcher, add the prefetcher name in the configuration file. 
Add the prefetcher to the configuration file like this.
```
{
    "L2C": {
        "prefetcher": "mypref"
    }
}
```
Prefetcher used in this project:
1. L1I: "FNLMMA"
2. L1D: "L1D_ipcp"
3. L2C: "L2C_ipcp"



```
$ ./config.sh <configuration file>
$ make
```

# Download QMM traces

You can download the QMM traces used in the project from here: 


# Run simulation

Execute the binary directly.
```
$ bin/champsim --warmup_instructions 30000000 --simulation_instructions 30000000 ~/path/to/traces/600.perlbench_s-210B.champsimtrace.xz
```


Execute using following script to run for all traces.
```
$ ./773TraceRun.sh
```

The number of warmup and simulation instructions given will be the number of instructions retired. Note that the statistics printed at the end of the simulation include only the simulation phase.



