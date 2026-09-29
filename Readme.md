Hello! Welcome to the code repository made for studying the long-range correlations and flow signatures in small systems.

The codebase is structured in the following way:

main/
├── run.sh                 # This shell script is tasked with running the event generation, post-processing, and plotting programs 
├── plots/                 # The finalised plots are stored here
├── data/                  # The necessary root files for the study 
├── generation/
│   └── gen.cc/         # Contains event generation code with:
|                          # PYTHIA tunes (default, ropes, shoving) + event loop:
│                          #   - cuts, PID applied 
│                          #   - HM/LM classification 
│                          #   - same-event Δφ pairing (3 subevent combos) 
│                          #   - fills correlation histograms directly (no mixed-event step for now)
├── post/                  # code made for calculations of Δφ projections, template fit (Eq. 2), V2Δ extraction, sqrt combination → v2(pT), and the plotting scripts    
├── Makefile               # auto-detects ROOT via root-config; requires $PYTHIA8 env var to be set      
└── README                 # You are here*


USER MANUAL: (You may make changes as you deem necessary)

Step 1: Execute the run.sh script using 'chmod +x run.sh', in the first run, then './run.sh' later, this line will automatically compile the codebase and begin the program execution

