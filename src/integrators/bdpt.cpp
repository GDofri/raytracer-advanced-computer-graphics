

#include <nori/camera.h>
#include <nori/rfilter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <nori/block.h>
#include <nori/timer.h>
#include <nori/bitmap.h>
#include <nori/gui.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#include <thread>
#include <utility>
#include "nori/bsdf.h"
#include "nori/emitter.h"
#include <cmath>
#include <filesystem>
#include <cstdlib>
#include <array>
NORI_NAMESPACE_BEGIN


class ThreadSafeVector {

public:
    ThreadSafeVector( size_t size ) : data(size) {
        for( size_t i = 0; i < size; i++ )
        {
            data.emplace_back(0.0f, 0);
        }
    }
    void add(const size_t &index, const float &value)
    {
        std::lock_guard<std::mutex> lock(m);
        data[index].first += value;
        data[index].second++;
    }
    void add(const std::vector<std::pair<float, int>> &values)
    {
        std::lock_guard<std::mutex> lock(m);
        for( size_t i = 0; i < values.size(); i++ )
        {
            data[i].first += values[i].first;
            data[i].second += values[i].second;
        }
    }
    std::vector<std::pair<float, int>> get()
    {
        std::lock_guard<std::mutex> lock(m);
        return data;
    }
    std::vector<std::pair<float, int>> data;
    std::mutex m;
};

/** Transport mode enum.
 */

enum class TransportMode {
    Radiance,
    Importance
};

/// BDPT utility function (From PBRTv3)
/// wo and wi can be swapped without affecting the outcome.
/// \param its The intersection
/// \param woWorld The Outgoing direction (in world coordinates, pointing away from its)
/// \param wiWorld The Incoming direction (in world coordinates, pointing away from its)
/// \param mode TransportMode

float correctShadingNormal(const Intersection &its, const Vector3f &woWorld,
                           const Vector3f &wiWorld, TransportMode mode) {

    if ( TransportMode::Importance == mode ) {
        float num   = std::abs(wiWorld.dot(its.shFrame.n)) * std::abs(woWorld.dot(its.geoFrame.n));
        float denom = std::abs(wiWorld.dot(its.geoFrame.n)) * std::abs(woWorld.dot(its.shFrame.n));
        if (denom == 0){
            return 0;
        }
        return num / denom;
    } else
        return 1;
}

struct Vertex {
    enum class Type {
        ECamera,
        ELight,
        ESurface
    };
public:
    Vertex() = default;

    /// Constructor for a vertex on the camera
    /// \param its Intersection data structure
    /// \param wi Incident direction in world coordinates
    /// \param alpha Path throughput
    /// \param isCamera True if the vertex is on the camera.
    Vertex(Intersection its, Vector3f wiWorld, Color3f alpha, Type type = Type::ESurface) :
        its(std::move(its)),
        wiWorld(std::move(wiWorld)),
        alpha(std::move(alpha)),
        pdfFwd(-1.f),
        pdfRev(-1.f),
        type(type)
    {};

    /// Constructor for a vertex on the path
    /// \param its Intersection data structure
    /// \param wi Incident direction in world coordinates
    /// \param alpha Path throughput
    /// \param pdfFwd Forward probability density
    /// \param isCamera True if the vertex is on the camera.
    Vertex(Intersection its, Vector3f wiWorld, Color3f alpha, float pdfFwd, Type type = Type::ESurface) :
        its(std::move(its)),
        wiWorld(std::move(wiWorld)),
        alpha(std::move(alpha)),
        pdfFwd(pdfFwd),
        pdfRev(-1.f),
        type(type)
    {};

    Color3f f(const Vector3f& p, TransportMode mode) const
    {
        if( Type::ECamera == type ){
            return 1;}
        Vector3f toP = p - its.p;
        if(toP.squaredNorm() == 0) return 0.f;
        Vector3f woWorld = toP.normalized();
        EMeasure measure = its.mesh->getBSDF()->isDiffuse() ? EMeasure::ESolidAngle : EMeasure::EDiscrete;
        BSDFQueryRecord bRec(its.toLocal(wiWorld), its.toLocal(woWorld), measure);
        Color3f f = its.mesh->getBSDF()->eval(bRec) * correctShadingNormal(its,
                                                                           woWorld,
                                                                           wiWorld,
                                                                           mode);
        return f;
    }

    // From PBRTv3 book:
    // The Vertex::Pdf() method returns the probability per unit area of the
    // sampling technique associated with a given vertex. Given a preceding
    // vertex prev, it evaluates the density for sampling the vertex next for
    // rays leaving the vertex *this. The prev argument may be equal to
    // nullptr for path endpoints (i.e., cameras or light sources), which
    // have no predecessor. Light sources require some extra care and are
    // handled separately via the PdfLight() method that will be discussed shortly.
    float pdf(const Scene &scene, const Vertex *prev,
              const Vertex &next) const {
        if( Type::ELight == type ) return PdfLight(next);
        float pdf;
        // Compute directions to preceding and next vertex
        Vector3f woWorld = next.its.p - its.p;
        if (woWorld.squaredNorm() == 0) return 0;
        woWorld.normalize();
        Vector3f wiWorld;

        if (prev) {
            // *this is a surface
            wiWorld = prev->its.p - its.p;
            if (wiWorld.squaredNorm() == 0) return 0;
            wiWorld.normalize();
            pdf = its.mesh->getBSDF()->pdf({its.toLocal(wiWorld), its.toLocal(woWorld), EMeasure::ESolidAngle});
        }
        else{
            // *this is a camera
            float unused;
            if( Type::ECamera != type ) { throw NoriException("pdf() called on non-camera vertex without prev"); }
            scene.getCamera()->pdf(Ray3f(its.p, woWorld), unused, pdf);
        }

        // Return probability per unit area at vertex _next_
        return convertDensity(pdf, next);
    }

    float PdfLight(const Vertex &v) const {

        Vector3f woWorld = ( v.its.p - its.p ).normalized();
        Emitter *emitter = its.mesh->getEmitter();
        EmitterQueryRecord eqr;
        eqr.woWorld = woWorld;
        eqr.normal = its.shFrame.n;
        float pdfDir = emitter->pdfDir(eqr);
        float cos = std::abs(v.its.shFrame.n.dot(woWorld));
        float geoTerm = cos/(its.p - v.its.p).squaredNorm();
        return pdfDir * geoTerm;
    }

    bool isConnectable() const {
        return Vertex::Type::ECamera == type
            || Vertex::Type::ELight == type
            || its.mesh->getBSDF()->isDiffuse();
    }

    bool isDelta() const {
        return ( !deltaOverride && ( Type::ECamera != type ) && !its.mesh->getBSDF()->isDiffuse() );
    }
    Intersection its;
    Vector3f wiWorld;
    Color3f alpha;

    // Pdfs in terms of area densities.
    float pdfFwd;
    float pdfRev;
    Type type;
    bool deltaOverride = false;

    // ConvertDensity funciton from PBRTv3
    float convertDensity(float pdf, const Vertex &next) const{

        Vector3f woWorld = next.its.p - its.p;
        float lengthSquared = woWorld.squaredNorm();
        woWorld.normalize();
        if( lengthSquared == 0) { return 0; }
        pdf *= std::abs(next.its.shFrame.n.dot(-woWorld));
        pdf /= lengthSquared;
        return pdf;
    }
};

// ScopedAssignment class from PBRTv3
template <typename Type>
class ScopedAssignment {
public:
    // ScopedAssignment Public Methods
    ScopedAssignment(Type *target = nullptr, Type value = Type())
            : target(target) {
        if (target) {
            backup = *target;
            *target = value;
        }
    }
    ~ScopedAssignment() {
        if (target) *target = backup;
    }
    ScopedAssignment(const ScopedAssignment &) = delete;
    ScopedAssignment &operator=(const ScopedAssignment &) = delete;
    ScopedAssignment &operator=(ScopedAssignment &&other) {
        if (target) *target = backup;
        target = other.target;
        backup = other.backup;
        other.target = nullptr;
        return *this;
    }

private:
    Type *target, backup;
};



class BDPTIntegrator : public Integrator {


    void printPath( const std::vector<Vertex> &path ) const
    {
        int i = 0;
        for( const Vertex &v : path )
        {
            std::cout << "Vertex " << i << ": ""Pos: " << v.its.p << " Wi: " << v.wiWorld << " Alpha: " << v.alpha << " pdf: " << v.pdfFwd << std::endl;
            // Flush
            std::cout << std::flush;
            i++;
        }
    }
    int generateLightSubpath( const Scene *scene, Sampler *sampler, std::vector<Vertex> &path, const Point2f& pixelSample, int maxDepth ) const
    {
        auto [emitter, emitterPdf] = scene->sampleEmitter(sampler->next1D());
        EmitterQueryRecord eqr;
        Ray3f ray;
        emitter->getEmitter()->sampleRay(eqr, *sampler, ray);

        // v1 (one indexed)
        // In veach : beta1 = W_e^(0)(y_0) / P_A(y_0)
        //          : p1 = P_A(y_0)
        Color3f beta = eqr.radiance; // As done in pbrt, no pdf in beta in first one/ (emitterPdf * eqr.pdfLoc);
        Intersection its;
        its.p = ray.o;
        its.shFrame = Frame(eqr.normal);
        its.mesh = emitter;
        path.emplace_back(its, its.shFrame.toLocal(-ray.d), beta, eqr.pdfLoc * emitterPdf, Vertex::Type::ELight); // As in pbrt, exclude emitterpdf except for next beta

        beta *= ray.d.dot(eqr.normal) / (emitterPdf * eqr.pdfDir * eqr.pdfLoc);

        return randomWalk( scene, sampler, path, maxDepth, eqr.pdfDir, beta, ray, TransportMode::Importance);
    }
    int generateCameraSubpath( const Scene *scene, Sampler *sampler, std::vector<Vertex> &path, const Point2f& pixelSample, int maxDepth ) const
    {
        Ray3f ray;
        // Sample a ray from the camera. Using -1,-1 as the aperture sample as it is not used in the camera implementation
        Point2f apertureSample = Point2f(-1.0f, -1.0f);
        Color3f beta = scene->getCamera()->sampleRay(ray, pixelSample, apertureSample);
        Intersection its;
        its.p = ray.o;
        its.shFrame = Frame(scene->getCamera()->getCameraNormal());

        // v1 (one indexed)
        // In veach : beta1 = W_e^(0)(z_0) / P_A(z_0) = 1/1 = 1
        //          : p1 = P_A(z_0) = 1
        // Here W_e^(0)(z_0) is the importance from a point on the camera...
        //
        path.emplace_back(its, Vector3f(-1.0f), beta,  1.0f, Vertex::Type::ECamera);

        float pdf = 1.0f;

        // In Veach: beta2 = f_s(z_-1 -> z_0 -> z_1) <-> z_1) * p1
        // Here f_s(z_-1 -> z_0 -> z_1) is W_e^(1)(z_0). It is the importance of the light from the camera to the
        // first vertex or as Veach calls it the directional component of the importance.
        // P_ortho(z_0 -> z_1) is the probability of sam* beta1 / P_ortho(z_0 -> z_1)
        //        //        //         : p2 = P_ortho(z_0 -> z_1) * G(z_0 pling the direction from the camera to the first vertex. Which
        float unused;
        scene->getCamera()->pdf(ray, unused, pdf );
        return randomWalk( scene, sampler, path, maxDepth, pdf, beta, ray, TransportMode::Radiance);
    }

    int randomWalk(const Scene *scene, Sampler *sampler, std::vector<Vertex> &path, const int maxDepth, const float pdf, Color3f beta, Ray3f ray, TransportMode mode ) const
    {
        float pdfFwd = pdf;
        float pdfRev = 0.f;
        int depth = 1;
        while( true )
        {
            Intersection its;
            if( !scene->rayIntersect( ray, its ) )
            {
                break;
            }
            Vector3f wi = its.toLocal(-ray.d);

            Vertex &prev = path[depth - 1];
            Vertex vertex = Vertex(its, -ray.d, beta);
            // pdfFwd converted into area density
            vertex.pdfFwd = prev.convertDensity(pdfFwd, vertex);
            path.push_back(vertex);

            if (++depth >= maxDepth) break;
            BSDFQueryRecord bRec(wi);
            bRec.measure = EMeasure::ESolidAngle;

            // in veach this is f_s(y_i-3 -> y_i-2 -> y_1)/P(y_i-2 <-> y_i-1)
            Color3f f = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            Vector3f wo = bRec.wo;
            pdfFwd = its.mesh->getBSDF()->isDiffuse() ? its.mesh->getBSDF()->pdf(bRec) : 1.0f;
            if (f.isZero() || pdfFwd == 0.f) break;

            // f includes both multiplication by the cosine term and division by the pdf
            beta *= f;

//            // Russian roulette
//            Removed as it brakes the assumption that pdfRev == pdfFwd
//            if (depth > 3) {
//                float rrProb = std::min(beta.maxCoeff(), 0.99f);
//                if (sampler->next1D() > rrProb) break;
//                beta /= rrProb;
//                pdfFwd *= rrProb;
//            }

            // Assuming that the pdf is the same for the reverse path
            pdfRev = pdfFwd;
            if( !its.mesh->getBSDF()->isDiffuse() ) {
                pdfFwd = 0.f;
                pdfRev = 0.f;
            }

            beta *= correctShadingNormal(its, its.toWorld(bRec.wo), its.toWorld(bRec.wi), mode);
            ray = Ray3f(its.p, its.toWorld(bRec.wo));
            prev.pdfRev = vertex.convertDensity(pdfRev, prev);
        }

        return depth;
    }
    Color3f connect( const Scene *scene, const std::vector<Vertex> &cameraPath, const std::vector<Vertex> &lightPath, int s, int t, Sampler *sampler, Point2f &lightRasterLoc, Vertex &sampled ) const
    {
        if( s == 0 )
        {
            // Check if last vertex on camera path is on a light source
            Vertex vertexT = cameraPath[t-1];
            if( vertexT.its.mesh->isEmitter() && vertexT.its.shFrame.n.dot(vertexT.wiWorld) > 0)
            {
                return vertexT.alpha * vertexT.its.mesh->getEmitter()->getRadiance();
            }
            return Color3f(0.0f);
        }
        if( t == 0 )
        {
            throw NoriException("t == 0 not implemented");
        }
        if( t==1 )
        {
            // Directly sample a camera instead of connecting to a vertex on the camera path
            const Vertex &vertexS = lightPath[s-1];
            if(!vertexS.isConnectable()) return 0.0f;
            Point3f cameraLoc = scene->getCamera()->getWorldPosition();
            float pdf;
            Color3f imp = scene->getCamera()->sampleWi(vertexS.its, lightRasterLoc, pdf);
            if( imp.isZero() || pdf <= 0.f ) return 0.0f;

            Vector3f vecToCamera = cameraLoc - vertexS.its.p;
            Intersection camIts;
            camIts.p = cameraLoc;
            camIts.shFrame = Frame(scene->getCamera()->getCameraNormal());
            sampled = Vertex(camIts, vecToCamera.normalized(), imp/pdf, 1.0f, Vertex::Type::ECamera);

            float cosI = std::abs(vecToCamera.normalized().dot(vertexS.its.shFrame.n));

            // This f does not include foreshortening factor.
            Color3f f = vertexS.f(cameraLoc, TransportMode::Importance);

            Color3f L = vertexS.alpha * f * sampled.alpha * cosI;
            if( L.isZero() ) return {0.0f};
            Ray3f shadowRay = Ray3f(vertexS.its.p, vecToCamera.normalized(), Epsilon, vecToCamera.norm() - Epsilon);
            if( scene->rayIntersect(shadowRay) )
            {
                return {0.0f};
            }
            return L;
        }
        if( s == 1 )
        {
            const Vertex &vertexT = cameraPath[t-1];
            if(!vertexT.isConnectable()) return 0.0f;

            // Directly sample a light source instead of connecting to a vertex on the light path
            auto [emitter, emitterPdf] = scene->sampleEmitter(sampler->next1D());
            EmitterQueryRecord eqr;
            emitter->getEmitter()->sample( eqr, *sampler);
            float lightLocPdf = eqr.pdfLoc;
            // Visability test
            Vector3f vecToLight = eqr.loc - vertexT.its.p;

            Intersection its;
            its.p = eqr.loc;
            its.mesh = emitter;
            its.shFrame = Frame(eqr.normal);
            sampled = Vertex( its, vecToLight.normalized(), eqr.radiance/(emitterPdf * eqr.pdfLoc), 0, Vertex::Type::ELight);
            sampled.pdfFwd = lightLocPdf*emitterPdf;

            float cosO = std::max(0.f,  eqr.normal.dot(-vecToLight.normalized()));
            float cosI = std::abs(vertexT.its.shFrame.n.dot(vecToLight.normalized()));
            if( cosI <= 0 || cosO <= 0 )
            {
                return Color3f(0.0f);
            }

            float g = cosO / vecToLight.squaredNorm();
            // This vertex.f does not include foreshortening factor so the cosI is needed.
            Color3f L = vertexT.alpha * vertexT.f( eqr.loc, TransportMode::Radiance) * sampled.alpha * g * cosI;

            Ray3f shadowRay = Ray3f(vertexT.its.p, vecToLight.normalized(), Epsilon, vecToLight.norm() - Epsilon);
            if ( L.isZero() || scene->rayIntersect(shadowRay))
            {
                return Color3f(0.0f);
            }
            return L;
        }

        // General case
        const Vertex &vertexS = lightPath[s-1];
        const Vertex &vertexT = cameraPath[t-1];
        if( !vertexS.isConnectable() || !vertexT.isConnectable() ) return 0.0f;
//        if( !vertexS.its.mesh->getBSDF()->isDiffuse() && !vertexT.its.mesh->getBSDF()->isDiffuse() )
//        {
//            return Color3f(0.0f);
//        }
        Color3f L = vertexS.alpha * vertexS.f(vertexT.its.p, TransportMode::Importance) * vertexT.f(vertexS.its.p, TransportMode::Radiance) * vertexT.alpha;
        if( L.isZero() ) return {0.0f};
        // Geometric term
        Vector3f inter = vertexT.its.p - vertexS.its.p;
        Ray3f shadowRay = Ray3f(vertexS.its.p, inter.normalized(), Epsilon, inter.norm() - Epsilon);
        if( scene->rayIntersect(shadowRay) )
        {
            return {0.0f};
        }
        Vector3f d = inter.normalized();
        float g = std::abs(vertexS.its.shFrame.n.dot(d))
                * std::abs(vertexT.its.shFrame.n.dot(-d))
                / inter.squaredNorm();
        return L * g;
    }

    size_t strategyIndex(int s, int t) const
    {
        int n = s+t;
        return s + (n + 1) * (n - 2) / 2;
    }

    // This is implemented in the fashion as done in pbrtv3
    float calculateMISWeight( std::vector<Vertex> &cameraPath, std::vector<Vertex> &lightPath, int s, int t, const Vertex &sampled, const Scene &scene ) const
    {
        // There are no alternative strategies for
        // cases of s = 2 or t = 2 (i.e. paths have single edge)
        if(s + t == 2) return 1;

        Vertex *vertexS = s > 0 ? &lightPath[s - 1] : nullptr;
        Vertex *vertexT = t > 0 ? &cameraPath[t - 1] : nullptr;
        Vertex *vertexSPrev = s > 1 ? &lightPath[s - 2] : nullptr;
        Vertex *vertexTPrev = t > 1 ? &cameraPath[t - 2] : nullptr;

        ScopedAssignment<Vertex> assignSampled;
        if(s == 1) { assignSampled = ScopedAssignment<Vertex>(vertexS, sampled); }
        if(t == 1) { assignSampled = ScopedAssignment<Vertex>(vertexT, sampled); }

        // Update how the last camera path vertex could have been sampled from the next vertex
        // as the next vertex has now changed from when the path was generated.
        ScopedAssignment<float> assignTPdfRev;
        if( s == 0) {
            float emitterPdf = scene.emitterSamplePdf();
            float pdfPos = vertexT->its.mesh->getEmitter()->pdfPos({});
            assignTPdfRev = {&vertexT->pdfRev,
                                                    emitterPdf * pdfPos};
        } else {
            assignTPdfRev = {&vertexT->pdfRev,
                                                     vertexS->pdf(scene, vertexSPrev, *vertexT)};
        }

        ScopedAssignment<float> assignTPrevPdfRev;
        if(vertexTPrev) {
            assignTPrevPdfRev = {&vertexTPrev->pdfRev,
                                                      s == 0 ? vertexT->PdfLight(*vertexTPrev) :
                                                      vertexT->pdf(scene, vertexS, *vertexTPrev)};
        }

        ScopedAssignment<float> assignSPdfRev;
        if(vertexS) {
            assignSPdfRev = {&vertexS->pdfRev,
                                                  vertexT->pdf(scene, vertexTPrev, *vertexS)};
        }

        ScopedAssignment<float> assignSPrevPdfRev;
        if(vertexSPrev) {
            assignSPrevPdfRev = {&vertexSPrev->pdfRev,
                                 vertexS->pdf(scene, vertexT, *vertexSPrev)};
        }

        // helper function to deal with delta bsdfs
        auto remap0 = [](float pdf){ return pdf == 0.f ? 1.f : pdf; };

        float riTotal = 0.f;
        float ri = 1.f;
        for (int i = t - 1; i > 0; --i) {
            ri *= remap0(cameraPath[i].pdfRev) / remap0(cameraPath[i].pdfFwd);
            if (!cameraPath[i].isDelta() && !cameraPath[i - 1].isDelta()) {
                riTotal += ri;
            }
        }

        ri = 1.f;
        for (int i = s - 1; i >= 0; --i) {
            ri *= remap0(lightPath[i].pdfRev) / remap0(lightPath[i].pdfFwd);
            if (!lightPath[i].isDelta() && ( (i==0) || !lightPath[i-1].isDelta() ) ){
                riTotal += ri;
            }
        }

        float weight = 1.f / (1.f + riTotal);
        return weight;
    }

    void renderBlock(const Scene *scene, Sampler *sampler, ImageBlock &block, ImageBlock &lightImage,
                     std::vector<std::shared_ptr<ImageBlock>> strategyBlocks, std::vector<std::pair<float, int>> &strategyWeights ) const {

        Point2i offset = block.getOffset();
        Vector2i size  = block.getSize();

        /* Clear the block contents */
        block.clear();

        int maxDepth = 40;
        std::vector<Vertex> cameraPath;
        cameraPath.reserve(maxDepth);
        std::vector<Vertex> lightPath;
        lightPath.reserve(maxDepth);

        /* For each pixel and pixel sample sample */
        for (int y=0; y<size.y(); ++y) {
            for (int x=0; x<size.x(); ++x) {
                for (uint32_t i=0; i<sampler->getSampleCount(); ++i) {

                    Point2f pixelSample = Point2f((float) (x + offset.x()), (float) (y + offset.y())) + sampler->next2D();
                    cameraPath.clear();
                    lightPath.clear();
                    int cameraPathLength = generateCameraSubpath(scene, sampler, cameraPath, pixelSample, maxDepth);
                    int lightPathLength = generateLightSubpath(scene, sampler, lightPath, pixelSample, maxDepth);

                    Color3f value = 0.0f;
                    for(int t = 1; t <= cameraPathLength; t++) {
                        for(int s = 0; s <= lightPathLength; s++) {
                            int depth = t + s -2;
                            if (t == 0 || (s == 1 && t == 1) || depth < 0 ||
                                depth > maxDepth)
                                continue;
                            Point2f lightRasterLoc;
                            Vertex sampled;
                            Color3f light = connect( scene, cameraPath, lightPath, s, t, sampler, lightRasterLoc, sampled );
                            float weight = 1;
                            if(!light.isZero()) {
                                weight = calculateMISWeight(
                                        cameraPath, lightPath, s, t, sampled, *scene);
                            }
                            light *= weight;

                            size_t strategy = strategyIndex(s, t);
                            if( !strategyBlocks.empty()  && strategy < strategyBlocks.size() )
                            {
                                strategyWeights[strategy].first += weight;
                                strategyWeights[strategy].second++;
                                if( t != 1){
                                    strategyBlocks[strategy]->put(pixelSample, light);
                                }
                                else {
                                    strategyBlocks[strategy]->splat(lightRasterLoc, light);
                                }
                            }

                            if( t != 1) {
                                value += light;
                            }
                            else if (!light.isZero()) {
                                lightImage.splat(lightRasterLoc, light);
                            }
                        }
                    }
                    /* Store in the image block */
                    block.put(pixelSample, value);
                }
            }
        }
    }

    void render(const Scene *scene, const std::string &fileName, int threadCount, bool showGUI) const {


        bool visualizeStrategies = true;
        showGUI = true;
        const Camera *camera = scene->getCamera();
        Vector2i outputSize = camera->getOutputSize();
        preprocess(scene);

        /* Create a block generator (i.e. a work scheduler) */
        BlockGenerator blockGenerator(outputSize, NORI_BLOCK_SIZE);

        /* Allocate memory for the entire output image and clear it */
        ImageBlock result(outputSize, camera->getReconstructionFilter());
        ImageBlock lightImage(outputSize, camera->getReconstructionFilter());
        result.clear();
        lightImage.clear();

        const int strategyFloors = 5;
        const int visualizationCount = (strategyFloors + 1) * (strategyFloors + 2) / 2 - 1;
        std::array<int, strategyFloors> lightImageStrategies = {};
        for( int i = 1; i < strategyFloors; i++)
        {
            lightImageStrategies[i] = (i+1)*(i+4)/2 - 1;
        }
        std::array<std::string, visualizationCount> strategyNames;
        ThreadSafeVector strategyWeights(visualizationCount);

        for( int floor = 1, stratIdx = 0; floor <= strategyFloors; floor++){
            for(int i = 0; i < floor + 1; ++i, ++stratIdx){
                std::string s = std::to_string(i);
                std::string t = std::to_string(floor+1-i);
                std::string d = std::to_string(floor-1);
                strategyNames[stratIdx] = "d" + d + "_s" + s + "_t" + t;
            }
        }

        std::vector<std::unique_ptr<ImageBlock>> strategyImages;
        if(visualizeStrategies)
        {
            for(int i = 0; i < visualizationCount; i++)
            {
                strategyImages.emplace_back(new ImageBlock(outputSize, camera->getReconstructionFilter()));
                strategyImages[i]->clear();
            }
        }

        /* Create a window that visualizes the partially rendered result */
        NoriScreen *screen = nullptr;
        if (showGUI) {
            nanogui::init();
            screen = new NoriScreen(result);
        }

        /* Do the following in parallel and asynchronously */
        std::thread render_thread([&] {
            tbb::task_scheduler_init init(threadCount);

            cout << "Rendering .. " << std::endl;
            cout.flush();
            Timer timer;

            tbb::blocked_range<int> range(0, blockGenerator.getBlockCount());

            auto map = [&](const tbb::blocked_range<int> &range) {
                /* Allocate memory for a small image block to be rendered
                   by the current thread */
                ImageBlock block(Vector2i(NORI_BLOCK_SIZE),
                                 camera->getReconstructionFilter());
                ImageBlock subLightImage(outputSize,
                                         camera->getReconstructionFilter());
                std::vector<std::shared_ptr<ImageBlock>> strategyBlocks;
                std::vector<std::pair<float, int>> strategyWeightsForThread;
                strategyWeightsForThread.reserve(visualizationCount);


                if(visualizeStrategies)
                {
                    for(int i = 0; i < visualizationCount; i++)
                    {
                        if( std::find(lightImageStrategies.begin(), lightImageStrategies.end(), i) ) {
                            // If strategy is a light image, make the block the whole output size.
                            strategyBlocks.emplace_back(new ImageBlock(outputSize, camera->getReconstructionFilter()));
                        }
                        else {
                            strategyBlocks.emplace_back(new ImageBlock(Vector2i(NORI_BLOCK_SIZE), camera->getReconstructionFilter()));
                        }
                        strategyBlocks[i]->clear();

                        strategyWeightsForThread.emplace_back(0,0);
                    }
                }

                /* Create a clone of the sampler for the current thread */
                std::unique_ptr<Sampler> sampler(scene->getSampler()->clone());

                for (int i=range.begin(); i<range.end(); ++i) {
                    /* Request an image block from the block generator */
                    blockGenerator.next(block);

                    /* Inform the sampler about the block to be rendered */
                    sampler->prepare(block);
                    if( visualizeStrategies )
                    {
                        for( std::shared_ptr<ImageBlock> &strategyBlock : strategyBlocks )
                        {
                            if( !std::find(lightImageStrategies.begin(), lightImageStrategies.end(), i) )
                            {
                                strategyBlock->setOffset(block.getOffset());
                                strategyBlock->setSize(block.getSize());
                            }
                        }
                    }
                    /* Render all contained pixels */
                    renderBlock(scene, sampler.get(), block, subLightImage, strategyBlocks, strategyWeightsForThread);

                    /* The image block has been processed. Now add it to
                       the "big" block that represents the entire image */
                    result.put(block);
                    lightImage.put(subLightImage);

                    strategyWeights.add(strategyWeightsForThread);

                    if(visualizeStrategies)
                    {
                        for( int i = 0; i < visualizationCount; i++)
                        {
                            strategyImages[i]->put(*strategyBlocks[i]);
                        }
                    }
                }
            };

            /// Default: parallel rendering
            tbb::parallel_for(range, map);

            /// (equivalent to the following single-threaded call)
//             map(range);

            cout << "done. (took " << timer.elapsedString() << ")" << endl;
        });

        /* Enter the application main loop */
        if (showGUI)
            nanogui::mainloop(50.f);

        /* Shut down the user interface */
        render_thread.join();

        if (showGUI) {
            delete screen;
            nanogui::shutdown();
        }

        ImageBlock combined(outputSize, camera->getReconstructionFilter());
        combined.clear();
        combined.put(result);
        combined.put(lightImage);
        /* Now turn the rendered image block into
           a properly normalized bitmap */
        float invSampleCount = 1.f / scene->getSampler()->getSampleCount();
        std::unique_ptr<Bitmap> bitmap(result.toBitmap(1.f));
        std::unique_ptr<Bitmap> lightBitmap(lightImage.toBitmap(invSampleCount));
        std::unique_ptr<Bitmap> combinedBitmap(combined.toBitmap(invSampleCount));

        auto createDirectory = [](const std::string path) {
            if( !std::filesystem::exists(path) ) {
                if( !std::filesystem::create_directories(path) ) {
                    std::cerr << "Failed to create directory " << path << std::endl;
                }
            }
        };

        /* Determine the filename of the output bitmap */
        size_t lastslash = fileName.find_last_of("/");
        size_t lastdot = fileName.find_last_of(".");
        std::string name = fileName.substr(lastslash + 1, lastdot - lastslash - 1);
        name += "_" + scene->getIntegrator()->getName() +
                      "_" + std::to_string(scene->getSampler()->getSampleCount());
        std::string directory = fileName.substr(0, lastslash + 1) + "results/" + name;
        std::string directoryExr = directory + "/exr";
        std::string directoryPng = directory + "/png";
        createDirectory(directoryExr);
        createDirectory(directoryPng);

        /* Save using the OpenEXR format */
        auto saveExr = [&directoryExr, &name] (Bitmap &image, const std::string &nameEnd) {
            image.saveEXR(directoryExr + "/" + name + nameEnd);
        };
        auto savePng = [&directoryPng, &name] (Bitmap &image, const std::string &nameEnd) {
            image.savePNG(directoryPng + "/" + name + nameEnd);
        };

        saveExr(*bitmap, "");
        saveExr(*lightBitmap, "_light");
        saveExr(*combinedBitmap, "_combined");
        savePng(*combinedBitmap, "_combined");

        if(visualizeStrategies)
        {
            for(int i = 0; i < visualizationCount; i++)
            {
                std::unique_ptr<Bitmap> strategyBitmap(strategyImages[i]->toBitmap(invSampleCount));
                saveExr(*strategyBitmap, "_" + strategyNames[i]);
                savePng(*strategyBitmap, "_" + strategyNames[i]);
            }
        }

//        auto weights = strategyWeights.get();
//        for( int i = 0; i < visualizationCount; i++)
//        {
//            float avgWeight = weights[i].first / weights[i].second;
//            cout << "Strategy " << strategyNames[i] << ": avg weight = " << avgWeight
//                 << " (" << weights[i].second << " samples)" << endl;
//        }

    }
    public:
    BDPTIntegrator(const PropertyList &props) {
    }
    std::string toString() const {
        return "BDPTIntegrator[]";
    }

    std::string getName() const {
        return "BDPT";
    }
};



    NORI_REGISTER_CLASS(BDPTIntegrator, "bdpt");
NORI_NAMESPACE_END