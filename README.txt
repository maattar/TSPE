Replication Codes for "Technology and Survival in Preindustrial England"

M. A. Attar

Note: Please refer to the comments within the *.m files for details.

--------------------------------------------------------------------------------------------

The first command to run is
>> run_Benchmark_Decennial(0,0,0,0,0)

This computes and saves the benchmark identifications of A_t and s_t

Then, the following commands should be executed:
>> dtests
>> irfs
>> mcycle

--------------------------------------------------------------------------------------------

For robustness, the identification code should be executed for different options:
>> run_Benchmark_Decennial(1,0,0,0,0)
>> run_Benchmark_Decennial(0,1,0,0,0) 
>> ...

Then, for identification with 25-year data, the following command should be executed:
>> run_Clark_Quadranscentennial

Finally, the following command generates the robustness results:
>> robustness

--------------------------------------------------------------------------------------------

For results using Broadberry et al. (2015) data, the following command should be executed:
>> run_Broadberry_Decennial
