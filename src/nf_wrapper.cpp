/// @file nf_wrapper.cxx
/// @brief C-compatible wrappers around the nanoflann KD-tree library.
///
/// This file provides a minimal C API for building a KD-tree over
/// column-major point data and performing k-nearest-neighbor queries.
/// The interface is intended for use from C or other languages via FFI.

#include <cassert>
#include <nanoflann.hpp>

#if defined _WIN32 || defined __CYGWIN__
  #define NF_EXPORT __declspec(dllexport)
#else
  #define NF_EXPORT __attribute__((visibility("default")))
#endif

/// @brief Adaptor for separate X and Y coordinate arrays (2D specialized).
///
/// @tparam T Scalar type of the coordinates (usually double).
template<typename T>
struct XYAdaptor {

    const int n_;
    const T *x_, *y_;

    /// @brief Constructor for separate arrays.
    /// @param d Unused in this adaptor, but kept for signature compatibility with NFHandleImpl.
    /// @param n Number of points.
    /// @param x Pointer to X coordinates.
    /// @param y Pointer to Y coordinates.
    XYAdaptor(int d, int n, const T* x, const T* y) :
        n_(n), x_(x), y_(y) {}

    /// @brief Return the number of points in the dataset.
    inline size_t kdtree_get_point_count() const { return n_; }

    /// @brief Access a single coordinate value.
    inline T kdtree_get_pt(const size_t idx, const size_t dim) const {
        return (dim == 0) ? x_[idx] : y_[idx];
    }

    /// @brief Bounding-box query (not implemented).
    template <class BBOX>
    bool kdtree_get_bbox(BBOX & /* bb */) const { return false; }
};


/// @brief Adaptor for separate X, Y, and Z coordinate arrays (3D specialized).
///
/// @tparam T Scalar type of the coordinates (usually double).
template<typename T>
struct XYZAdaptor {

    const int n_;
    const T *x_, *y_, *z_;

    /// @brief Constructor for separate X, Y, Z arrays.
    /// @param d Unused, kept for signature compatibility with NFHandleImpl.
    /// @param n Number of points.
    /// @param x Pointer to X coordinates.
    /// @param y Pointer to Y coordinates.
    /// @param z Pointer to Z coordinates.
    XYZAdaptor(int d, int n, const T* x, const T* y, const T* z) :
        n_(n), x_(x), y_(y), z_(z) {}

    /// @brief Return the number of points in the dataset.
    inline size_t kdtree_get_point_count() const { return n_; }

    /// @brief Access a single coordinate value.
    /// Note: nanoflann internal loops for DIM=3 will optimize this switch.
    inline T kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0) return x_[idx];
        if (dim == 1) return y_[idx];
        return z_[idx];
    }

    /// @brief Bounding-box query (returning false defaults to internal calculation).
    template <class BBOX>
    bool kdtree_get_bbox(BBOX & /* bb */) const { return false; }
};


/// @brief Minimal nanoflann adaptor for column-major point data.
///
/// @tparam T Scalar type of the coordinates.
template<typename T>
struct SimpleAdaptor {

    const int d_, n_;
    const T* data_;

    // N.b. caller must guarantee lifetime of data
    SimpleAdaptor(int d, int n, T* data) :
        d_(d), n_(n), data_(data) {}

    /// @brief Return the number of points in the dataset.
    inline size_t kdtree_get_point_count() const { return n_; }

    /// @brief Access a single coordinate value.
    ///
    /// @param idx Point index.
    /// @param dim Coordinate index.
    inline T kdtree_get_pt(const size_t idx, const size_t dim) const {
        return data_[dim + idx*d_];
    }

    /// @brief Bounding-box query (not implemented).
    template <class BBOX>
    bool kdtree_get_bbox(BBOX & /* bb */) const { return false; }

};

// Hidden from the C header, only in your .cpp
struct NFHandle {
    virtual ~NFHandle() = default;
    virtual void knn_query(double* query_pt, int nn, int* index, double* dist) = 0;
};

template<typename T_Adaptor, int DIM = -1>
struct NFHandleImpl final : public NFHandle {
    using TreeType = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, T_Adaptor>,
        T_Adaptor, DIM, int>;

    T_Adaptor adaptor;
    TreeType tree;

    // Perfect forwarding of constructor args to the adaptor
    template<typename... Args>
    NFHandleImpl(int dims, Args&&... args)
        : adaptor(dims, std::forward<Args>(args)...),
          tree(dims, adaptor) {}

    void knn_query(double* query_pt, int nn, int* index, double* dist) override {
        // Direct call to knnSearch writes straight into the C-provided buffers
        [[maybe_unused]] size_t nfound =
            tree.knnSearch(query_pt, static_cast<size_t>(nn), index, dist);
        assert(nfound == static_cast<size_t>(nn));
    }
};

extern "C" {

/// @brief Create and initialize a KD-tree.
///
/// @param dims Number of dimensions.
/// @param n    Number of points.
/// @param data Pointer to point coordinates (column-major).
///
/// @return Opaque handle to the KD-tree, or nullptr on failure.
///
/// @note The input data must remain valid until @ref c_nf_free is called.
NF_EXPORT void *c_nf_init(int dims, int n, double *data) {
    if (dims <= 0 || n <= 0 || !data) return nullptr;

    NFHandle *ptr = nullptr;

    // Specialize for common dimensions to trigger compiler optimizations
    if (dims == 2) {
        ptr = new NFHandleImpl<SimpleAdaptor<double>, 2>(2, n, data);
    } else if (dims == 3) {
        ptr = new NFHandleImpl<SimpleAdaptor<double>, 3>(3, n, data);
    } else {
        // Fallback for any other dimension (runtime-defined)
        ptr = new NFHandleImpl<SimpleAdaptor<double>>(dims, n, data);
    }

    return ptr;
}


/// @brief Create and initialize a 2D KD-tree from separate X and Y arrays.
///
/// @note The input data must remain valid until @ref c_nf_free is called.
NF_EXPORT void* c_nf_init_xy(int n, double* x, double* y) {
    if (n <= 0 || !x || !y) return nullptr;

    // We pass '2' for both the template DIM and the runtime dims parameter
    return new NFHandleImpl<XYAdaptor<double>, 2>(2, n, x, y);
}


/// @brief Create and initialize a 3D KD-tree from separate X, Y, and Z arrays.
///
/// @note The input data must remain valid until @ref c_nf_free is called.
NF_EXPORT void* c_nf_init_xyz(int n, double* x, double* y, double* z) {
    if (n <= 0 || !x || !y || !z) return nullptr;

    // We pass '3' for both the template DIM and the runtime dims parameter.
    // The variadic constructor in NFHandleImpl forwards (n, x, y, z) to the XYZAdaptor.
    return new NFHandleImpl<XYZAdaptor<double>, 3>(3, n, x, y, z);
}


/// @brief Destroy a KD-tree and release associated resources.
///
/// @param ptr Opaque handle returned by @ref c_nf_init or @ref c_nf_init_xy.
NF_EXPORT void c_nf_free(void *ptr) {
    if (!ptr) return;

    // Cast to the base interface.
    // The virtual destructor ensures the correct subclass (~NFHandleImpl) is called.
    NFHandle *h = static_cast<NFHandle *>(ptr);
    delete h;
}


/// @brief Perform a n-nearest-neighbor query.
///
/// @param ptr  Opaque KD-tree handle.
/// @param p    Query point (array of length @c dims).
/// @param nn   Number of nearest neighbors to find.
/// @param idx  Output array of length @c nn with point indices.
/// @param dist Output array of length @c nn with squared distances.
///
/// @pre The KD-tree must have been created with @ref c_nf_init.
/// @pre Output arrays must be preallocated by the caller.
NF_EXPORT
void c_nf_n_nearest(void *ptr, double p[], int nn, int idx[], double dist[]) {

    // Cast to the interface instead of the concrete nf_handle
    auto* instance = static_cast<NFHandle*>(ptr);

    // This call dispatches to the correct NFHandleImpl<Adaptor, DIM>::knn_query
    instance->knn_query(p, nn, idx, dist);
}

} // extern "C"








