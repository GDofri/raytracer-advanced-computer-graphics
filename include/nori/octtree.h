//
// Created by dofri on 29.2.2024.
//

#ifndef NORI_OCTTREE_H
#define NORI_OCTTREE_H

#include <memory>
#include "common.h"
#include "bbox.h"
#include "mesh.h"
#include <unordered_set>
#include <array>
#include <optional>
NORI_NAMESPACE_BEGIN

struct Node{
    nori::BoundingBox3f bbox;
    std::array<Node*, 8> children;
    std::unordered_set<uint32_t> triangles;
    bool isLeaf;
};

class Octtree
{
    public:


        Octtree(const nori::BoundingBox3f& bbox, const nori::Mesh &mesh);
        std::unordered_set<uint32_t> getPossibleIntersections(const Ray3f &ray);
        Node *getRoot();

        /** getTreeStatistics
         *  Computes some characteristics of the tree.
         * @return A tuple of three integers. The first one is the number of leaf nodes,
         * the second is the number of interior nodes,
         * the third the collective number of triangles in all the leaf nodes (duplicates included).
         */
        std::tuple<uint32_t, uint32_t, uint32_t> getTreeStatistics();

        /**
         *
         * @param ray
         * @return
         */
        std::pair< uint32_t, std::optional<Intersection> > findClosestIntersection(const Ray3f &ray, bool isShadowRay = false);

private:

    /**
     * Constructs the octtree.
     */
    Node *build(const BoundingBox3f& bbox, const std::unordered_set<uint32_t> &triangleIndices, uint32_t);
    std::unordered_set<uint32_t> getPossibleIntersections(const Ray3f &ray, Node *node, uint32_t depth);


    Node *m_root;
    const nori::MatrixXu &m_indices;
    const Mesh &m_mesh;

    static const uint32_t MAXDEPTH = 8;
    bool triangleIntersectsBoundingBox(uint32_t triangleIdx, const nori::BoundingBox3f &bbox);


    void traverseForStatistics(Node *node, uint32_t &interior, uint32_t &leaf, uint32_t &triangles);

    std::pair< uint32_t, std::optional<Intersection>> findClosestIntersectionFromSet(const Ray3f &ray, const std::unordered_set<uint32_t> &triangles, bool isShadowRay);

    std::pair< uint32_t, std::optional<Intersection> > traverseForClosestIntersection(const Ray3f &ray, Node *node, uint32_t depth, bool isShadowRay);
};


#endif //NORI_OCTTREE_H
NORI_NAMESPACE_END