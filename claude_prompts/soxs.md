# Develop the NTT/SOXS spectrograph

## Goals

We wish to develop the NTT/SOXS spectrograph for PypeIt.  

## Claude

### Skills

Consider using the skills in PypeIt-development-suite/.claude/skills/

## Context

Here are the context documents:
- Webpage for the instrument: https://www.eso.org/public/teles-instr/lasilla/ntt/soxs/
- SOXS pipeine: /home/xavier/Projects/PypeIt/Others/soxspipe

The PypeIt pipeline

## Coding

Here are guidelines for coding: 

- Use Python where possible
- Add inline comments to explain the effort
- Reuse existing code when possible
- Place import statements at the top of the file.
- Include a description of inputs/outputs in the doc string of all methods
- Add "Created by JXP and Claude" to any new methods.

## Testing

If you need to test the code: 

- use the files in PypeIt-development-suite/pypeitdev/ntt_soxs/
- run on the "pypeit14b" conda environment.

## Prompts

### Spectrograph Class

1. Read the context documents to understand the instrument and the code.  Your first task is to create a new spectrograph class in PypeIt/pypeit/spectrographs/ntt_soxs.py.  Model it after the other spectrograph classes in that folder.  Use the files in PypeIt-development-suite/pypeitdev/ntt_soxs/ as a guide.  Before doing so, ask me questions about the instrument and the code in Q&A below.  Use Fable if you can.  Log your work

## Q&A

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs