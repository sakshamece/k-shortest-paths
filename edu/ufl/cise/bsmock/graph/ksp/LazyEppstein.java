package edu.ufl.cise.bsmock.graph.ksp;

import edu.ufl.cise.bsmock.graph.*;
import edu.ufl.cise.bsmock.graph.util.*;
import java.util.*;

/**
 * Lazy version of Eppstein's algorithm (by Jimenez and Marzal) for computing the K shortest paths
 * between two nodes in a graph.
 *
 * Created by Brandon Smock on October 5, 2015.
 * Last updated by Brandon Smock on December 24, 2015.
 */
public final class LazyEppstein extends Eppstein implements KSPAlgorithm {

    public LazyEppstein() {}

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
        HashMap<String, EppsteinArrayHeap> arrayHeaps = new HashMap<>(graph.numNodes());

        EppsteinHeap hg = new EppsteinHeap(new Edge(sourceLabel, sourceLabel, 0, 0.0));

        ArrayList<Path> ksp = new ArrayList<>();
        PriorityQueue<EppsteinPath> pathPQ = new PriorityQueue<>();

        pathPQ.add(new EppsteinPath(hg, -1, tree.getNodes().get(sourceLabel).getDist()));

        for (int i = 0; i < K && pathPQ.size() > 0; i++) {
            EppsteinPath kpathImplicit = pathPQ.poll();
            Path kpath = kpathImplicit.explicitPath(ksp, tree);
            ksp.add(kpath);
            addExplicitChildrenToQueue(kpathImplicit, ksp, pathPQ);
            if (!outrootHeaps.containsKey(kpathImplicit.getHeap().getSidetrack().getToNode())) {
                buildHeap(kpathImplicit.getHeap().getSidetrack().getToNode(), graph, sidetrackEdgeCostMap, nodeHeaps, edgeHeaps, outrootHeaps, arrayHeaps, tree);
            }
            addCrossEdgeChildToQueue(outrootHeaps, kpathImplicit, i, ksp, pathPQ);
        }

        return ksp;
    }

    private static void buildHeap(String nodeLabel, Graph graph, HashMap<String, Double> sidetrackEdgeCostMap, HashMap<String, EppsteinHeap> nodeHeaps, HashMap<String, EppsteinHeap> edgeHeaps, HashMap<String, EppsteinHeap> outrootHeaps, HashMap<String, EppsteinArrayHeap> arrayHeaps, ShortestPathTree tree) {
        computeOutHeap(nodeLabel, graph, sidetrackEdgeCostMap, nodeHeaps, edgeHeaps);

        EppsteinArrayHeap currentArrayHeap;
        if (nodeLabel.equals(tree.getRoot())) {
            currentArrayHeap = new EppsteinArrayHeap();
        } else {
            if (!outrootHeaps.containsKey(tree.getParentOf(nodeLabel))) {
                buildHeap(tree.getParentOf(nodeLabel), graph, sidetrackEdgeCostMap, nodeHeaps, edgeHeaps, outrootHeaps, arrayHeaps, tree);
            }
            currentArrayHeap = arrayHeaps.get(tree.getParentOf(nodeLabel));
        }

        EppsteinHeap sidetrackHeap = nodeHeaps.get(nodeLabel);
        if (sidetrackHeap != null) {
            currentArrayHeap = currentArrayHeap.clone();
            currentArrayHeap.addOutroot(sidetrackHeap);
        }

        arrayHeaps.put(nodeLabel, currentArrayHeap);
        EppsteinHeap currentHeap = currentArrayHeap.toEppsteinHeap2();
        if (currentHeap != null) {
            outrootHeaps.put(nodeLabel, currentHeap);
        }
    }
}
