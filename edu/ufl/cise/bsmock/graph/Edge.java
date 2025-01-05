package edu.ufl.cise.bsmock.graph;

/**
 * The Edge class implements standard properties and methods for a weighted edge in a directed graph.
 *
 * Created by Brandon Smock on 6/19/15.
 */
public class Edge implements Cloneable {
    private String fromNode;
    private String toNode;
    private double weight;
    private double availableBW; // New field

    public Edge() {
        this.fromNode = null;
        this.toNode = null;
        this.weight = Double.MAX_VALUE;
        this.availableBW = 0.0; // Initialize new field
    }

    public Edge(String fromNode, String toNode, double weight, double availableBW) {
        this.fromNode = fromNode;
        this.toNode = toNode;
        this.weight = weight;
        this.availableBW = availableBW; // Initialize new field
    }

    public String getFromNode() {
        return fromNode;
    }

    public void setFromNode(String fromNode) {
        this.fromNode = fromNode;
    }

    public String getToNode() {
        return toNode;
    }

    public void setToNode(String toNode) {
        this.toNode = toNode;
    }

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }

    public double getAvailableBW() {
        return availableBW;
    }

    public void setAvailableBW(double availableBW) {
        this.availableBW = availableBW;
    }

    public Edge clone() {
        return new Edge(fromNode, toNode, weight, availableBW);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("(");
        sb.append(fromNode);
        sb.append(",");
        sb.append(toNode);
        sb.append("){");
        sb.append(weight);
        sb.append(",");
        sb.append(availableBW); // Append new field
        sb.append("}");

        return sb.toString();
    }

    public boolean equals(Edge edge2) {
        if (hasSameEndpoints(edge2) && weight == edge2.getWeight() && availableBW == edge2.getAvailableBW())
            return true;

        return false;
    }

    public boolean hasSameEndpoints(Edge edge2) {
        if (fromNode.equals(edge2.getFromNode()) && toNode.equals(edge2.getToNode()))
            return true;

        return false;
    }
}
