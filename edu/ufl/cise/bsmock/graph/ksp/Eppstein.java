package edu.ufl.cise.bsmock.graph.ksp;

import edu.ufl.cise.bsmock.graph.*;
import edu.ufl.cise.bsmock.graph.util.*;
import java.util.*;

/**
 * Eppstein's algorithm for computing the K shortest paths between two nodes in a graph.
 *
 * Created by Brandon Smock on October 5, 2015.
 * Last updated by Brandon Smock on December 24, 2015.
 */
public class Eppstein implements KSPAlgorithm {

    public boolean isLoopless() {
        return false;
    }

    public Eppstein() {}

    public List<Path> ksp(Graph graph, String sourceLabel, String targetLabel, int K) {
        ShortestPathTree tree;
        try {
            tree = Dijkstra.shortestPathTree(graph.transpose(), targetLabel);
        } catch (Exception e) {
            tree = new ShortestPathTree(targetLabel);
        }

        HashMap<String,Double> sidetrackEdgeCostMap = computeSidetrackEdgeCosts(graph, tree);

        HashMap<String,EppsteinHeap> nodeHeaps = new HashMap<>(graph.numNodes());
        HashMap<String,EppsteinHeap> edgeHeaps = new HashMap<>(graph.numEdges());
        HashMap<String,EppsteinHeap> outrootHeaps = new HashMap<>();

        for (String nodeLabel : graph.getNodes().keySet()) {
            computeOutHeap(nodeLabel, graph, sidetrackEdgeCostMap, nodeHeaps, edgeHeaps);
        }

        Graph reversedSPT = new Graph();
        for (DijkstraNode node: tree.getNodes().values()) {
            reversedSPT.addEdge(node.getParent(), node.getLabel(), graph.getNode(node.getLabel()).getNeighbors().get(node.getParent()));
        }

        EppsteinArrayHeap rootArrayHeap = new EppsteinArrayHeap();
        recursiveOutrootHeaps(targetLabel, rootArrayHeap, nodeHeaps, outrootHeaps, reversedSPT);

        EppsteinHeap hg = new EppsteinHeap(new Edge(sourceLabel, sourceLabel, 0, 0.0));

        ArrayList<Path> ksp = new ArrayList<>();
        PriorityQueue<EppsteinPath> pathPQ = new PriorityQueue<>();

        pathPQ.add(new EppsteinPath(hg, -1, tree.getNodes().get(sourceLabel).getDist()));

        for (int i = 0; i < K && pathPQ.size() > 0; i++) {
            EppsteinPath kpathImplicit = pathPQ.poll();
            Path kpath = kpathImplicit.explicitPath(ksp, tree);
            ksp.add(kpath);
            addExplicitChildrenToQueue(kpathImplicit, ksp, pathPQ);
            addCrossEdgeChildToQueue(outrootHeaps, kpathImplicit, i, ksp, pathPQ);
        }

        return ksp;
    }

    protected static HashMap<String,Double> computeSidetrackEdgeCosts(Graph graph, ShortestPathTree tree) {
        HashMap<String, Double> sidetrackEdgeCostMap = new HashMap<>();
        List<Edge> edgeList = graph.getEdgeList();
        for (Edge edge : edgeList) {
            String tp = tree.getParentOf(edge.getFromNode());
            if (tp == null || !tp.equals(edge.getToNode())) {
                double sidetrackEdgeCost = edge.getWeight() + tree.getNodes().get(edge.getToNode()).getDist() - tree.getNodes().get(edge.getFromNode()).getDist();
                sidetrackEdgeCostMap.put(edge.getFromNode() + "," + edge.getToNode(), sidetrackEdgeCost);
            }
        }
        return sidetrackEdgeCostMap;
    }

    protected static void computeOutHeap(String nodeLabel, Graph graph, HashMap<String,Double> sidetrackEdgeCostMap, HashMap<String,EppsteinHeap> nodeHeaps, HashMap<String,EppsteinHeap> edgeHeaps) {
        Node node = graph.getNode(nodeLabel);
        ArrayList<Edge> sidetrackEdges = new ArrayList<>();
        Edge bestSidetrack = null;
        double minSidetrackCost = Double.MAX_VALUE;
        for (String neighbor : node.getAdjacencyList()) {
            String edgeLabel = nodeLabel + "," + neighbor;
            if (sidetrackEdgeCostMap.containsKey(edgeLabel)) {
                double sidetrackEdgeCost = sidetrackEdgeCostMap.get(edgeLabel);
                if (sidetrackEdgeCost < minSidetrackCost) {
                    if (bestSidetrack != null) {
                        sidetrackEdges.add(bestSidetrack);
                    }
                    bestSidetrack = new Edge(nodeLabel, neighbor, node.getNeighbors().get(neighbor), 0.0); // Initialize with availableBW
                    minSidetrackCost = sidetrackEdgeCost;
                } else {
                    sidetrackEdges.add(new Edge(nodeLabel, neighbor, node.getNeighbors().get(neighbor), 0.0)); // Initialize with availableBW
                }
            }
        }

        if (bestSidetrack != null) {
            EppsteinHeap bestSidetrackHeap = new EppsteinHeap(bestSidetrack, sidetrackEdgeCostMap.get(bestSidetrack.getFromNode() + "," + bestSidetrack.getToNode()));
            EppsteinArrayHeap arrayHeap = new EppsteinArrayHeap();
            if (sidetrackEdges.size() > 0) {
                bestSidetrackHeap.setNumOtherSidetracks(bestSidetrackHeap.getNumOtherSidetracks() + 1);
                for (Edge edge : sidetrackEdges) {
                    EppsteinHeap sidetrackHeap = new EppsteinHeap(edge, sidetrackEdgeCostMap.get(edge.getFromNode() + "," + edge.getToNode()));
                    edgeHeaps.put(edge.getFromNode() + "," + edge.getToNode(), sidetrackHeap);
                    arrayHeap.add(sidetrackHeap);
                }
                bestSidetrackHeap.addChild(arrayHeap.toEppsteinHeap());
            }
            nodeHeaps.put(nodeLabel, bestSidetrackHeap);
            edgeHeaps.put(bestSidetrack.getFromNode() + "," + bestSidetrack.getToNode(), bestSidetrackHeap);
        }
    }

    protected static void addExplicitChildrenToQueue(EppsteinPath kpathImplicit, ArrayList<Path> ksp, PriorityQueue<EppsteinPath> pathPQ) {
        double kpathCost = kpathImplicit.getCost();
        for (EppsteinHeap childHeap : kpathImplicit.getHeap().getChildren()) {
            int prefPath = kpathImplicit.getPrefPath();
            Double candidateCost = ksp.get(prefPath).getTotalCost() + childHeap.getSidetrackCost();
            EppsteinPath candidate = new EppsteinPath(childHeap, prefPath, candidateCost);
            pathPQ.add(candidate);
        }
    }

    protected static void addCrossEdgeChildToQueue(HashMap<String,EppsteinHeap> outrootHeaps, EppsteinPath kpathImplicit, int prefPath, ArrayList<Path> ksp, PriorityQueue<EppsteinPath> pathPQ) {
        if (outrootHeaps.containsKey(kpathImplicit.getHeap().getSidetrack().getToNode())) {
            EppsteinHeap childHeap = outrootHeaps.get(kpathImplicit.getHeap().getSidetrack().getToNode());
            Double candidateCost = ksp.get(prefPath).getTotalCost() + childHeap.getSidetrackCost();
            EppsteinPath candidate = new EppsteinPath(childHeap, prefPath, candidateCost);
            pathPQ.add(candidate);
        }
    }

    protected static void recursiveOutrootHeaps(String nodeLabel, EppsteinArrayHeap currentArrayHeap, HashMap<String,EppsteinHeap> nodeHeaps, HashMap<String,EppsteinHeap> outrootHeaps, Graph reversedSPT) {
        EppsteinHeap sidetrackHeap = nodeHeaps.get(nodeLabel);
        if (sidetrackHeap != null) {
            currentArrayHeap = currentArrayHeap.clone();
            currentArrayHeap.addOutroot(sidetrackHeap);
        }
        EppsteinHeap currentHeap = currentArrayHeap.toEppsteinHeap2();
        if (currentHeap != null) {
            outrootHeaps.put(nodeLabel, currentHeap);
        }
        for (String neighbor : reversedSPT.getNode(nodeLabel).getNeighbors().keySet()) {
            recursiveOutrootHeaps(neighbor, currentArrayHeap, nodeHeaps, outrootHeaps, reversedSPT);
        }
    }
}

class EppsteinHeap {
    private Edge sidetrack;
    private double sidetrackCost = 0.0;
    private ArrayList<EppsteinHeap> children;
    private int numOtherSidetracks = 0;

    public EppsteinHeap(Edge sidetrack) {
        this.sidetrack = sidetrack;
        this.children = new ArrayList<>();
    }

    public EppsteinHeap(Edge sidetrack, Double sidetrackCost) {
        this.sidetrack = sidetrack;
        this.sidetrackCost = sidetrackCost;
        this.children = new ArrayList<>();
    }

    public EppsteinHeap(Edge sidetrack, double sidetrackCost, ArrayList<EppsteinHeap> children, int numOtherSidetracks) {
        this.sidetrack = sidetrack;
        this.sidetrackCost = sidetrackCost;
        this.children = children;
        this.numOtherSidetracks = numOtherSidetracks;
    }

    public Edge getSidetrack() {
        return sidetrack;
    }

    public void setSidetrack(Edge sidetrack) {
        this.sidetrack = sidetrack;
    }

    public double getSidetrackCost() {
        return sidetrackCost;
    }

    public void setSidetrackCost(double sidetrackCost) {
        this.sidetrackCost = sidetrackCost;
    }

    public ArrayList<EppsteinHeap> getChildren() {
        return children;
    }

    public void setChildren(ArrayList<EppsteinHeap> children) {
        this.children = children;
    }

    public void addChild(EppsteinHeap child) {
        this.children.add(child);
    }

    public int getNumOtherSidetracks() {
        return numOtherSidetracks;
    }

    public void setNumOtherSidetracks(int numOtherSidetracks) {
        this.numOtherSidetracks = numOtherSidetracks;
    }

    public EppsteinHeap clone() {
        ArrayList<EppsteinHeap> children_clone = new ArrayList<>(children.size());
        for (EppsteinHeap eh: children) {
            children_clone.add(eh);
        }
        return new EppsteinHeap(sidetrack, sidetrackCost, children_clone, numOtherSidetracks);
    }
}

class EppsteinArrayHeap {
    private ArrayList<EppsteinHeap> arrayHeap;

    public EppsteinArrayHeap() {
        arrayHeap = new ArrayList<>(0);
    }

    public ArrayList<EppsteinHeap> getArrayHeap() {
        return arrayHeap;
    }

    public void setArrayHeap(ArrayList<EppsteinHeap> arrayHeap) {
        this.arrayHeap = arrayHeap;
    }

    public int getParentIndex(int i) {
        return (i-1)/2;
    }

    public void add(EppsteinHeap h) {
        arrayHeap.add(h);
        bubbleUp(arrayHeap.size()-1);
    }

    public void addOutroot(EppsteinHeap h) {
        int current = arrayHeap.size();
        while (current > 0) {
            int parent = getParentIndex(current);
            EppsteinHeap newHeap = arrayHeap.get(parent).clone();
            arrayHeap.set(parent, newHeap);
            current = parent;
        }
        arrayHeap.add(h);
        bubbleUp(arrayHeap.size() - 1);
    }

    private void bubbleUp(int current) {
        if (current == 0) return;
        int parent = getParentIndex(current);
        if (arrayHeap.get(current).getSidetrackCost() >= arrayHeap.get(parent).getSidetrackCost()) return;
        EppsteinHeap temp = arrayHeap.get(current);
        arrayHeap.set(current, arrayHeap.get(parent));
        arrayHeap.set(parent, temp);
        bubbleUp(parent);
    }

    public EppsteinHeap toEppsteinHeap() {
        if (arrayHeap.size() == 0) return null;
        EppsteinHeap eh = arrayHeap.get(0);
        for (int i = 1; i < arrayHeap.size(); i++) {
            EppsteinHeap h = arrayHeap.get(i);
            arrayHeap.get(getParentIndex(i)).addChild(h);
        }
        return eh;
    }

    public EppsteinHeap toEppsteinHeap2() {
        int current = arrayHeap.size() - 1;
        if (current == -1) return null;
        while (current >= 0) {
            EppsteinHeap childHeap = arrayHeap.get(current);
            while (childHeap.getChildren().size() > childHeap.getNumOtherSidetracks()) {
                childHeap.getChildren().remove(childHeap.getChildren().size() - 1);
            }
            int child1 = current * 2 + 1;
            int child2 = current * 2 + 2;
            if (child1 < arrayHeap.size()) arrayHeap.get(current).addChild(arrayHeap.get(child1));
            if (child2 < arrayHeap.size()) arrayHeap.get(current).addChild(arrayHeap.get(child2));
            if (current > 0) current = getParentIndex(current);
            else current = -1;
        }
        return arrayHeap.get(0);
    }

    public EppsteinArrayHeap clone() {
        EppsteinArrayHeap clonedArrayHeap = new EppsteinArrayHeap();
        for (EppsteinHeap heap: arrayHeap) {
            clonedArrayHeap.add(heap);
        }
        return clonedArrayHeap;
    }
}

class EppsteinPath implements Comparable<EppsteinPath> {
    EppsteinHeap heap;
    int prefPath;
    Double cost;

    public EppsteinPath(EppsteinHeap heap, int prefPath, Double cost) {
        this.heap = heap;
        this.prefPath = prefPath;
        this.cost = cost;
    }

    public int getPrefPath() {
        return prefPath;
    }

    public void setPrefPath(int prefPath) {
        this.prefPath = prefPath;
    }

    public EppsteinHeap getHeap() {
        return heap;
    }

    public void setHeap(EppsteinHeap heap) {
        this.heap = heap;
    }

    public Double getCost() {
        return cost;
    }

    public void setCost(Double cost) {
        this.cost = cost;
    }

    public Path explicitPath(List<Path> ksp, ShortestPathTree tree) {
        Path explicitPath = new Path();
        if (prefPath >= 0) {
            Path explicitPrefPath = ksp.get(prefPath);
            LinkedList<Edge> edges = explicitPrefPath.getEdges();
            int lastEdgeNum = -1;
            Edge heapSidetrack = heap.getSidetrack();
            for (int i = edges.size() - 1; i >= 0; i--) {
                Edge currentEdge = edges.get(i);
                if (currentEdge.getToNode().equals(heapSidetrack.getFromNode())) {
                    lastEdgeNum = i;
                    break;
                }
            }
            explicitPath = new Path();
            for (int i = 0; i <= lastEdgeNum; i++) {
                explicitPath.add(edges.get(i));
            }
            explicitPath.add(heap.getSidetrack());
        }
        String current = heap.getSidetrack().getToNode();
        while (!current.equals(tree.getRoot())) {
            String next = tree.getParentOf(current);
            Double edgeWeight = tree.getNodes().get(current).getDist() - tree.getNodes().get(next).getDist();
            explicitPath.add(new Edge(current, next, edgeWeight, 0.0)); // Initialize with availableBW
            current = next;
        }
        return explicitPath;
    }

    public int compareTo(EppsteinPath comparedNode) {
        double cost1 = this.cost;
        double cost2 = comparedNode.getCost();
        if (cost1 == cost2) return 0;
        if (cost1 > cost2) return 1;
        return -1;
    }
}
