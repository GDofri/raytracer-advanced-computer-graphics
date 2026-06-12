//
// Created by dofri on 29.2.2024.
//

#include <utility>
#include <numeric>

#include "nori/octtree.h"
#include "nori/mesh.h"

NORI_NAMESPACE_BEGIN

Octtree::Octtree(const BoundingBox3f& bbox, const Mesh &mesh):
    m_indices(mesh.getIndices()),
    m_mesh(mesh)
{
    std::unordered_set< uint32_t > indexes(m_indices.cols());
    for( uint32_t i = 0; i < m_indices.cols(); i++)
    {
        indexes.insert(i);
    }
    m_root = build(bbox, indexes, 0);
}

Node *Octtree::build(const BoundingBox3f& bbox, const std::unordered_set<uint32_t> &triangleIndices, uint32_t depth)
{
    uint32_t numberOfTriangles = triangleIndices.size();

    if( 0 == numberOfTriangles)
    {
        return nullptr;
    }

    Node* node = new Node;
    node->bbox = bbox;
    node->isLeaf = ( numberOfTriangles < 10 || MAXDEPTH == depth );


    if( node->isLeaf ){
        node->triangles = triangleIndices;
        return node;
    }

    // Create sub-bounding-boxes
    auto min = bbox.min;
    auto max = bbox.max;
    auto ctr = bbox.getCenter();

    std::vector< BoundingBox3f > bboxes
            {
                {{min[0], min[1], min[2]},{ctr[0], ctr[1], ctr[2]}},
                {{ctr[0], min[1], min[2]},{max[0], ctr[1], ctr[2]}},
                {{min[0], ctr[1], min[2]},{ctr[0], max[1], ctr[2]}},
                {{ctr[0], ctr[1], min[2]},{max[0], max[1], ctr[2]}},
                {{min[0], min[1], ctr[2]},{ctr[0], ctr[1], max[2]}},
                {{ctr[0], min[1], ctr[2]},{max[0], ctr[1], max[2]}},
                {{min[0], ctr[1], ctr[2]},{ctr[0], max[1], max[2]}},
                {{ctr[0], ctr[1], ctr[2]},{max[0], max[1], max[2]}},
            };
    std::vector< std::unordered_set<uint32_t> > trianglesInBbox(8);

    // Find what triangles go into each child node
    for( const auto triangleIdx : triangleIndices)
    {
        for( size_t i = 0; i < 8; ++i)
        {
            if(triangleIntersectsBoundingBox(triangleIdx, bboxes[i]))
            {
                trianglesInBbox[i].insert(triangleIdx);
            }
        }
    }

    // Create children
    for( size_t i = 0; i < 8; ++i)
    {
        node->children[i] = build(bboxes[i], trianglesInBbox[i], depth + 1);
    }
    return node;
}

Node *Octtree::getRoot()
{
    return m_root;
}

bool Octtree::triangleIntersectsBoundingBox( uint32_t triangleIdx, const BoundingBox3f &bbox )
{
    return bbox.overlaps(m_mesh.getBoundingBox(triangleIdx));
}


std::unordered_set<uint32_t> Octtree::getPossibleIntersections(const Ray3f &ray)
{
    return getPossibleIntersections(ray, m_root, 0);
}


std::unordered_set<uint32_t> Octtree::getPossibleIntersections(const Ray3f &ray, Node *node, uint32_t depth)
{

    if( !node->bbox.rayIntersect(ray) )
    {
        return {};
    }
    if( node->isLeaf)
    {
        auto res = findClosestIntersectionFromSet(ray, node->triangles, false);
        if(res.first != (uint32_t)-1) return {res.first};
        else return {};
//        return node->triangles;
    }

    std::unordered_set<uint32_t> resultSet = {};
    for( const auto &child : node->children)
    {
        if(child)
        {
            auto childIntersections = getPossibleIntersections(ray, child, depth + 1);
            resultSet.insert(childIntersections.begin(), childIntersections.end());
        }
    }
    return resultSet;
}

std::pair< uint32_t, std::optional<Intersection> > Octtree::findClosestIntersection(const Ray3f &ray, bool isShadowRay)
{
    return traverseForClosestIntersection(ray, m_root, 0, isShadowRay);
}


std::pair< uint32_t, std::optional<Intersection>>  Octtree::traverseForClosestIntersection(const Ray3f &ray, Node *node, uint32_t depth, bool isShadowRay)
{
    if( !node->bbox.rayIntersect(ray) )
    {
        return {0, std::nullopt};
    }
    if( node->isLeaf)
    {
        // Find closest triangle intersection
        return findClosestIntersectionFromSet(ray, node->triangles, isShadowRay);
    }

    std::array<std::pair<Node*, float>, 8> childBBoxTs = {};
    for(size_t i = 0; i < 8; i++)
    {
        float t, _;
        auto child = node->children[i];

        if( child && child->bbox.rayIntersect(ray) && child->bbox.rayIntersect(ray, t, _))
        {
            childBBoxTs[i] = {child, t};
        }
        else
        {
            childBBoxTs[i] = {child, std::numeric_limits<float>::max()};
        }
    }


    auto byDistanceToIntersection = [](const std::pair<Node*, float> &childA, const std::pair<Node*, float> &childB)
    {
        return  childA.second < childB.second;
    };

    std::sort(childBBoxTs.begin(), childBBoxTs.end(), byDistanceToIntersection);

    uint32_t closestIdx = -1;
    Intersection closestIntersection;
    closestIntersection.t = std::numeric_limits<float>::max();
    bool foundIntersectionInChildren = false;
    for( const auto &[child, childT] : childBBoxTs)
    {
        // Because of the sorting, we can stop early if we find a child that we don't collide with.
        if( childT == std::numeric_limits<float>::max() ) break;

        if(child && childT < closestIntersection.t)
        {
            auto [closestChildIdx, intersection] = traverseForClosestIntersection(ray, child, depth + 1, isShadowRay);
            if( intersection.has_value() && intersection->t < closestIntersection.t )
            {
                closestIntersection = *intersection;
                closestIdx = closestChildIdx;
                foundIntersectionInChildren = true;
                if(isShadowRay) break;
            }
        }
    }
    return {closestIdx, foundIntersectionInChildren ? std::make_optional(closestIntersection) : std::nullopt};

}


std::pair< uint32_t, std::optional<Intersection>> Octtree::findClosestIntersectionFromSet(const Ray3f &ray, const std::unordered_set<uint32_t> &triangles, bool isShadowRay)
{
    uint32_t minIdx = std::numeric_limits<uint32_t>::max();
    float minT = std::numeric_limits<float>::max();
    Intersection intersection;

    bool foundIntersection = false;
    for( const auto &idx : triangles)
    {
        float u, v, t;
        if (m_mesh.rayIntersect(idx, ray, u, v, t) && t < minT)
        {
            intersection.uv = Point2f(u, v);
            intersection.mesh = &m_mesh;
            intersection.t = t;
            minIdx = idx;
            minT = t;
            foundIntersection = true;
            if(isShadowRay) break;
        }
    }

    return  {minIdx, foundIntersection ? std::make_optional(intersection) : std::nullopt};
}

void Octtree::traverseForStatistics( Node* node, uint32_t &interior, uint32_t &leaf, uint32_t &triangles)
{
    // Node is a leaf node
    if( node->isLeaf )
    {
        leaf++;
        triangles += node->triangles.size();
        return;
    }
    // Node is an interior node
    interior++;
    // Traverse child nodes
    for( const auto child : node->children)
    {
        if( !child ) continue;
        traverseForStatistics(child, interior, leaf, triangles);
    }
}

std::tuple<uint32_t, uint32_t, uint32_t> Octtree::getTreeStatistics()
{
    uint32_t interior = 0, leaf = 0, triangles = 0;

    traverseForStatistics(m_root, interior, leaf, triangles);
    return  {interior, leaf, triangles};
}
NORI_NAMESPACE_END



















