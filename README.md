# Area Law Violations in Entanglement Measures for Disordered and Inhomogeneous Quantum Spin Chains

## Intro
This repo contains my Master's Thesis completed under the supervision of Dr. Paola Ruggiero at King's College London (KCL) in September 2022. 

This work has not been peer reviewed and should be treated as such. I am sharing it because I have been asked to by students also working on the Strong Disorder Renormalization Group (SDRG) approach to spin chains. 

Comments and feedback are welcome. 

If you have any questions, please let me know!

## How to Use
Clone the repo:

```bash
git clone https://github.com/josephcbradley/sdrgjl
```

Instantiate the environment:
    
```bash
chmod +x instantiate.sh 
./instantiate.sh
```

Run the tests:

```bash
chmod +x runtests.sh
./runtests.sh
```

If everything looks OK, you can run the whole project:

```bash
chmod +x full_run.sh
./full_run.sh
```

This will produce all the figures and compile the LaTeX.

Use ```part_run.sh``` to run a subset of the project. 

I have tested this on a Macbook Pro with an Intel chip in 2024 on Julia 1.10.4 and everything seems to run. LaTeX is the most likely thing to break IMO. 

Have fun!