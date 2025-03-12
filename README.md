# LinearGraph-NormalTree

In this repository we implement the ideas discussed in [this documet from MIT](https://web.mit.edu/2.151/www/Handouts/EqFormulation.pdf) where they discuss how to describe physical systems as graphs and then how to get the equations describing the system. <br> 

We will begin by showing an example <br>

<center>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="/assets/tests/test17/fig17.svg" >
  <img alt="Figure 17'" src="/assets/tests/test17/fig17_light.svg"  width="300" height="300">
</picture>
</center>

<br>

Which is equivalent to the following graph

<br>

<center>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="/assets/tests/test17/graph17.svg" >
  <img alt="Graph 17'" src="/assets/tests/test17/graph17_light.svg" width="300" height="300">
</picture>
</center>

<br>

Then we build the normal tree of the graph, following this procedure:
- Include all across sources
- Include all possible "A" types (those included become independent components)
- Include all possible "D" types
- Include all possible "T" types (those not included become independent components)
- Check that no through sources can be included <br>

That results in the following tree:

<br>

<center>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="/assets/tests/test17/tree17.svg" >
  <img alt="Normal Tree 17'" src="/assets/tests/test17/tree17_light.svg" width="300" height="300">
</picture>
</center>

<br>


Now we generate the state equations using the following procedure:
- Writes $B-S$ elemental equations for all non-source edges
- Writes $N-1-S_A$ continuity equations by getting all the child nodes from both source and target nodes and assuring the flux in between is constant
- Writes $B-N+1-S_T$ compatibility equations by inserting a link into the tree and assuring the path from source to target in the tree equals the link

Now we reduce to the minimal equations only involving the independent state variables and the sources using Gaussian elimination. <br>

The result is a system $\dot{x} = Ax + Bs$ where $x$ is a vector written in a file `var_index.csv`, $A$ is written in `var.csv`, $s$ is written in `sources_index.csv` and $B$ is written in `sources.csv`. <br>

We can check that in this case we get 
$$x=(F_{K_1}\quad v_m)\quad \mathrm{ and }\quad s = (\dot{F_s}\quad F_s)$$ <br>

If you want to run the algorithm you must create a new folder `foldername` in `/assets` and a file `foldername/edges.json` with the information about the graph. Then put `foldername` next to `./Main` in the `run.txt` file. After, you can run `bash run.txxt` in a command line to get the results.




