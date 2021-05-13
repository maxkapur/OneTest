# Paper outline

## Background
General school-choice literature has focused on bilateral matching markets with arbitrary preferences

Or single-sided matching markets

In the bilateral case, the most important theoretical work is Azevedo and Leshno. Show the connection between stable matches and market clearing eq. Insight has been underappreciated.

However, that paper offers a theoretical description without considering the conditions under which equilibrium might be reached or computable instances.

The one computable instance it does consider is with completely iid preferences on the part of the schools. then it provides expressions for comparative statics and a theoretical discussion about when negative incentives can arise.

I propose a continuous model that is both computational tractable, moderately realistic, and versatile enough to accommodate not just centralized markets that used a DA procedure but also dynamic markets in which schools freely admit and reject students.

## Preliminaries
Proof of Azevedo and Leshno's theorem and its implications re: dynamic markets.
Characterizations of equilibrium:
- variational inequality
- convex program

The market used in this paper
- consideration sets stack
- characterization of demand function
- characterization of equilibrium

Argue for the relative realisticness of this market, real-world examples, linear combinations. 

## Solution strategies
Tatonnement procedure and convergence proof. 

Direct solution when cutoff ordering is given: show that $p = [A^{-1} (q - \gamma)]^+$

Heuristics for determining optimal cutoff ordering; potentially yields a closed-form solution. 

## Optimization tasks
Situations in which you can use only a subset of the parameters.

Kinds of optimization tasks and how to do them.

## Comparative statics
Using the direct solution above, how to derive them.

Edge cases when two $p_c$ are equal; argue about interpretation

## Extensions
Propose my general form of the dynamic admissions market.

Its complexity; the number of potential preference lists and consideration sets.

Computable instances include mine and iid scores.

Admissions coalitions and clusters. 