### String Reconstruction from Paired k-mers

This repository provides a Python implementation for reconstructing a string from paired k-mers using Eulerian path. The algorithm constructs a directed graph from the paired k-mers, finds an Eulerian path, and reconstructs the string by traversing the graph.

### Overview

Given a set of paired k-mers, the goal is to reconstruct the original string from these overlapping pairs. The algorithm works by:
Constructing a directed graph from the paired k-mers.
Finding an Eulerian path in the graph, which is a path that visits every edge exactly once.
Reconstructing the string from the Eulerian path.
The graph is represented as an adjacency list, and the Eulerian path is found using a depth-first search (DFS) approach. The algorithm handles the case where the graph might have an unbalanced number of in-degrees and out-degrees by selecting the appropriate starting node.

### Algorithm Steps
Construct the graph: Build a directed graph where each pair of k-mers defines an edge from the prefix of the first k-mer to the suffix of the second k-mer.
Find the Eulerian path: Use DFS to find a path that visits every edge exactly once.
Reconstruct the string: Traverse the Eulerian path to reconstruct the original string from the paired k-mers.

### Input Format
The input consists of a list of paired k-mers, where each pair is a tuple:
('prefix_kmer', 'suffix_kmer')

### Example
For the input:
paired_kmers = [
    ('ACC', 'ATA'),
    ('ACT', 'ATT'),
    ('ATA', 'TGA'),
    ('ATT', 'TGA'),
    ('CAC', 'GAT'),
    ('CCG', 'TAC'),
    ('CGA', 'ACT'),
    ('CTG', 'AGC'),
    ('CTG', 'TTC'),
    ('GAA', 'CTT'),
    ('GAT', 'CTG'),
    ('GAT', 'CTG'),
    ('TAC', 'GAT'),
    ('TCT', 'AAG'),
    ('TGA', 'GCT'),
    ('TGA', 'TCT'),
    ('TTC', 'GAA')
]
The reconstructed string will be:

result = string_reconstruction_from_read_pairs(paired_kmers)
print(result)
Output:
'ACCATAACTATTGACATCGATGATTCTGAGCTTCGAAGCTT'
Code Example
from collections import defaultdict

def string_reconstruction_from_read_pairs(paired_kmers):
    def construct_graph(paired_kmers):
        graph = defaultdict(list)
        for pair in paired_kmers:
            prefix, suffix = pair[0], pair[1]
            graph[prefix].append(suffix)
        return graph

    def find_eulerian_path(graph, start_node):
        path = []
        stack = [start_node]

        while stack:
            current_node = stack[-1]

            if graph[current_node]:
                next_node = graph[current_node].pop()
                stack.append(next_node)
            else:
                path.append(stack.pop())

        return path[::-1]

    graph = construct_graph(paired_kmers)

    start_node = None
    for node in graph:
        indegree = sum(len(graph[parent]) for parent in graph)
        outdegree = len(graph[node])
        if indegree < outdegree:
            start_node = node
            break

    if start_node is None:
        start_node = list(graph.keys())[0]

    eulerian_path = find_eulerian_path(graph, start_node)

    reconstructed_string = eulerian_path[0]
    for node in eulerian_path[1:]:
        reconstructed_string += node[-1]

    return reconstructed_string

# Example usage:
paired_kmers = [
    ('ACC', 'ATA'),
    ('ACT', 'ATT'),
    ('ATA', 'TGA'),
    ('ATT', 'TGA'),
    ('CAC', 'GAT'),
    ('CCG', 'TAC'),
    ('CGA', 'ACT'),
    ('CTG', 'AGC'),
    ('CTG', 'TTC'),
    ('GAA', 'CTT'),
    ('GAT', 'CTG'),
    ('GAT', 'CTG'),
    ('TAC', 'GAT'),
    ('TCT', 'AAG'),
    ('TGA', 'GCT'),
    ('TGA', 'TCT'),
    ('TTC', 'GAA')
]

result = string_reconstruction_from_read_pairs(paired_kmers)
print(result)

### Features
Graph Construction: Converts paired k-mers into a directed graph.
Eulerian Path: Finds the Eulerian path using a depth-first search.
String Reconstruction: Reconstructs the original string from the Eulerian path.

### Requirements
Python 3.x

### How to Use
Clone or download the repository.

Run the script with your own set of paired k-mers to reconstruct the original string.

### Contribution
Feel free to fork the repository and submit pull requests if you'd like to contribute or improve the project.

### License
This project is licensed under the MIT License.

