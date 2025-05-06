# Single-Source Shortest Path  (SSSP)

This project implements a **parallel algorithm** to update **single-source shortest paths (SSSP)** in large-scale **dynamic graphs**, based on the research paper:  
 **"A Parallel Algorithm Template for Updating Single-Source Shortest Paths in Large-Scale Dynamic Networks"**  
 *Authors: Khanda, Srinivasan, Bhowmick, Norris, Das (IEEE TPDS, 2022)*

---

## Objective

To design and implement a **scalable parallel algorithm** using **OpenMP** that efficiently updates SSSP trees in dynamic graphs where edge insertions and deletions occur — without recomputing everything from scratch.

---

## Technologies & Tools

- **C++** – Language used for implementation  
- **OpenMP** – For shared-memory, intra-node parallelism  
- **Dynamic Scheduling** – For workload balancing during updates  
- **CSR Format** – For graph representation  
- **Real-world Datasets** – From Network Repository and synthetic R-MAT graphs  
- **GitHub** – For collaboration and progress tracking

---

## Project Flow

1. **Graph Input**: Read static graph in CSR or adjacency list format  
2. **Edge Updates**: Load batch of inserted/deleted edges  
3. **Subgraph Detection**: Mark affected vertices  
4. **Parallel Updates**: Relax distances using OpenMP threads  
5. **Asynchronous Rounds**: Optionally control level of synchronization  
6. **Final Output**: Updated SSSP tree and performance stats

---

## Key Concepts from the Research

- **Dynamic Graphs**: Graphs where edges are frequently inserted or deleted
- **SSSP Tree Representation**: Uses a rooted tree storing distance and parent per vertex
- **Affected Subgraph Identification**: Only vertices impacted by changes are processed
- **Iterative Updates Without Locks**: Convergence via repeated distance checks instead of synchronization
- **Parallelism**:
  - Intra-node using **OpenMP**
  - Dynamic scheduling for load balancing
- **Avoid Cycle Formation**: During insertions, maintain acyclic SSSP tree
- **Asynchronous Updates**: Reduces synchronization cost by processing neighbors multiple levels deep
- **Batch Processing**: Handles massive edge changes in manageable chunks for performance

---

## Parallelization Strategy

## 1. MPI (Inter-node Communication)

- Could be used for distributing graph partitions (using METIS) across different nodes.
- Each node can independently update its subgraph’s SSSP tree.
- After local updates, MPI_Gather or MPI_Allreduce could be used to sync global distances.

## 2. OpenMP (Intra-node Parallelism)

- Already used in shared-memory version (OpenMP).
- Enables parallel processing of:
- Edge deletions/insertions
- Affected vertex updates

Dynamic scheduling in OpenMP helps in load balancing during updates of uneven subtrees.

## 3. METIS (Graph Partitioning)

- Can pre-process large graphs to partition them into smaller subgraphs.
- Each subgraph can be assigned to a thread (OpenMP) or node (MPI).
- Reduces cross-node communication and optimizes memory locality.

---

![Image](https://github.com/user-attachments/assets/c0644cdc-1f3f-4795-89c9-8d480f60ce13)


![Image](https://github.com/user-attachments/assets/6311abb0-2fff-48ff-9918-f3eddc5ade11)


![Image](https://github.com/user-attachments/assets/8c97e62c-0fea-4c09-863b-c859c0c7facc)


![Image](https://github.com/user-attachments/assets/6ce47c2f-e3a0-468f-ad8e-cee77c8898a5)


![Image](https://github.com/user-attachments/assets/accb118e-1e3b-48e9-b671-059fd78bf23e)


![Image](https://github.com/user-attachments/assets/0fcb7ec6-caae-42e4-bb82-f789b6f78a89)


![Image](https://github.com/user-attachments/assets/66e15b41-7aaf-408e-99dd-7d80de0ec23e)


![Image](https://github.com/user-attachments/assets/f38dd282-0b41-40eb-8d33-6d4c7cfb95b0)


![Image](https://github.com/user-attachments/assets/f518de16-ea10-48ce-b30f-9f727b8a87d8)


![Image](https://github.com/user-attachments/assets/4aefca2c-2d5b-4433-93e5-b57625682360)

## 👨‍💻 Contributors

- Hammad Shabbir  
- Haider Zia 
- Muhammad Iqrash Qureshi
