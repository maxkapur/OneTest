# Paper outline

## Intro
The classical *school choice problem* is a one-to-many stable assignment problem. Students match with schools; try to make a stable match.

Recently, attention has shifted toward a nonatomic formulation of this problem, in which individual students are replaced with a distribution of students over the space of possible preference lists and scores. 

This enables the characterization of assignment policies, including stable assignments, via school cutoffs. 
Azevedo and Leshno were the first to make this suggestion. This is useful, because most admissions markets are not run by a central agency. Instead, they are dynamic: colleges can admit or reject students as they please. School cutoffs make more sense as a decision variable in a dynamic market, because this is actually something schools can control, than assignments.

The cutoff formulation also creates space for the possibility that schools may have goals other than filling their capacity. They may not even have a meaningful limit on capacity, as is the case for many online schools.

In this dynamic situation, that paper offers a theoretical description of the relationship between supply and demand and posits a few comparative statics. However, it doesn't concern itself with computation, neither of the location of the equilibrium nor of the statics themselves. 

The one computable instance it does consider is with completely iid preferences on the part of the schools. then it provides expressions for comparative statics and a theoretical discussion about when negative incentives can arise.

It is this computational question to which I turn. I propose a continuous model that is both computationally tractable, moderately realistic, and versatile enough to accommodate not just centralized markets that used a DA procedure but also dynamic markets in which schools freely admit and reject students.

This model enables us to actually compute useful comparative statics. 

## Preliminaries
Parameters eta and q

Definition of equilibrium and possible interpretations. 

List possible optimization tasks. 


## Model
Characterization of students: MNL choice

Characterization of schools: same preference order. Argue for realism via linear combinations and market segmentation.

Later theoretical results require each school to have a capacity; we will use this later.

Demand function: derivation from stack of consideration sets.
Appeal function.
Matrix expression for each of the above; note that it depends on ordering of p.

## Optimization tasks
Compute the demand, and show how to actually do the following opt tasks 

### Computing the equilibrium
Equilibrium exists and is unique by ....

Tatonnement procedure and convergence proof.

When sum of capacities < 1, observe that D = q.

Direct solution when cutoff ordering is given: show that $p = [A^{-1} (q - \gamma)]^+$ if given the optimal ordering. 

Heuristic for determining optimal cutoff ordering.

Validate the results using deferred acceptance. 

### Reverse optimization of student preferences
Compute gamma; show inductive procedure. Argue for the informational power of this gamma. 

## Comparative statics
When p, D, and gamma are available:
### Cutoff effects
On demand, appeal when quality is fixed
### Quality effects
On demand, appeal when cutoffs are fixed

### At equilibrium: Demand effects
On cutoff when quality is fixed

Effect of a change in total student population


## Extensions
Propose my general form of the dynamic admissions market.

Equivalent problems to equilibrium:
- variational inequality
- convex program

Its complexity; the number of potential preference lists and consideration sets.

Computable instances include mine and iid scores.

Admissions coalitions and clusters. 