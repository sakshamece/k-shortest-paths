package edu.ufl.cise.bsmock.graph;

/**
 * The Node class implements a node in a directed graph keyed on a label of type String, with adjacency lists for
 * representing edges.
 *
 * Created by brandonsmock on 5/31/15.
 */

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;

public class Node {
    protected String label;
    protected HashMap<String, Double> neighbors; // adjacency list, with HashMap for each edge weight
    protected HashMap<String, Double> availableBWs; // adjacency list, with HashMap for each edge available bandwidth

    public Node() {
        neighbors = new HashMap();
        availableBWs = new HashMap();
    }

    public Node(String label) {
        this.label = label;
        neighbors = new HashMap();
        availableBWs = new HashMap();
    }

    public String getLabel() {
        return label;
    }

    public void setLabel(String label) {
        this.label = label;
    }

    public HashMap<String, Double> getNeighbors() {
        return neighbors;
    }

    public void setNeighbors(HashMap<String, Double> neighbors) {
        this.neighbors = neighbors;
    }

    public HashMap<String, Double> getAvailableBWs() {
        return availableBWs;
    }

    public void setAvailableBWs(HashMap<String, Double> availableBWs) {
        this.availableBWs = availableBWs;
    }

    public void addEdge(String toNodeLabel, Double weight, Double availableBW) {
        neighbors.put(toNodeLabel, weight);
        availableBWs.put(toNodeLabel, availableBW);
    }

    public double removeEdge(String toNodeLabel) {
        if (neighbors.containsKey(toNodeLabel)) {
            double weight = neighbors.get(toNodeLabel);
            neighbors.remove(toNodeLabel);
            availableBWs.remove(toNodeLabel);
            return weight;
        }

        return Double.MAX_VALUE;
    }

    public Set<String> getAdjacencyList() {
        return neighbors.keySet();
    }

    public LinkedList<Edge> getEdges() {
        LinkedList<Edge> edges = new LinkedList<Edge>();
        for (String toNodeLabel : neighbors.keySet()) {
            edges.add(new Edge(label, toNodeLabel, neighbors.get(toNodeLabel), availableBWs.get(toNodeLabel)));
        }

        return edges;
    }

    public String toString() {
        StringBuilder nodeStringB = new StringBuilder();
        nodeStringB.append(label);
        nodeStringB.append(": {");
        Set<String> adjacencyList = this.getAdjacencyList();
        Iterator<String> alIt = adjacencyList.iterator();
        HashMap<String, Double> neighbors = this.getNeighbors();
        while (alIt.hasNext()) {
            String neighborLabel = alIt.next();
            nodeStringB.append(neighborLabel.toString());
            nodeStringB.append(": ");
            nodeStringB.append(neighbors.get(neighborLabel));
            if (alIt.hasNext())
                nodeStringB.append(", ");
        }
        nodeStringB.append("}");
        nodeStringB.append("\n");

        return nodeStringB.toString();
    }
}
