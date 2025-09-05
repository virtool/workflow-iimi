# workflow-iimi

**This workflow is experimental**.

An analysis workflow for detecting known OTUs (viruses) using read mapping and ML.

## Steps

1. **Build all-OTU index**. Build a Bowtie2 mapping index that includes all sequences in 
   the ML model bundle.
2. **Map all OTUs**. Map sample reads against all OTUs in the ML model bundle reference.
3. **Predict**. Use Iimi (pre-release) to predict viral presence based on mappings.
