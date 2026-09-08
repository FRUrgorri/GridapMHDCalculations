using MPI

MPI.Init() #This command creates the fork

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

#Define total iterations and range for each process
total_iter = 100
chunk = div(total_iter, size) #The integer part of the division
start_idx = rank * chunk + 1
end_idx = (rank == size - 1) ? total_iter : (rank + 1) * chunk

#Perform work for this process's chunk
local_results = []
for i in start_idx:end_idx
    push!(local_results, i)
end
 
#Gather results in rank 0 (optional)

# If results are isbit (numbers) Gather can be used
"""
all_results = MPI.Gather(last(local_results), comm; root=0)
if rank == 0
    println("All results in MPI rank 0: ", all_results)
end
"""

# For other objects (like arrays in this example)

if rank == 0
    all_results = []  #To allocate the memory
    push!(all_results,local_results) # The firs process is the first element
    
    for i in 1:(size-1) #For over the rest of the processes
        push!(all_results, MPI.recv(comm; source=i, tag=10)) #10 is a tag of the sending operation
    end
    
    println(all_results)
else
    MPI.send(local_results,comm;dest=0,tag=10)    
end

MPI.Finalize()

