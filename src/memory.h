/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef _XMIPP_MEMORY
#define _XMIPP_MEMORY

#include <cstddef>
#include <cstdlib>

#include "src/error.h"

inline std::size_t xmipp_index_count(int lower, int upper)
{
    if (upper < lower)
        return 0;
    return static_cast<std::size_t>(
        static_cast<long long>(upper) - static_cast<long long>(lower) + 1);
}

/** Indexed pointer: logical v[i] maps to data[i - lower] inside the allocation. */
template <class T>
class XmippIndexPtr
{
public:
    T* data;
    int lower;

    XmippIndexPtr() : data(NULL), lower(0) {}
    XmippIndexPtr(T* data_ptr, int lower_bound) : data(data_ptr), lower(lower_bound) {}

    T& operator[](int i) const
    {
        return data[static_cast<std::ptrdiff_t>(i) - lower];
    }

    explicit operator bool() const { return data != NULL; }
    bool operator==(std::nullptr_t) const { return data == NULL; }
    bool operator!=(std::nullptr_t) const { return data != NULL; }
};

/** Packed Numerical Recipes 1-based flat matrix view over 0-based storage. */
template <class T>
class XmippIndexFlat2D
{
public:
    T* data;
    int cols;

    XmippIndexFlat2D() : data(NULL), cols(0) {}
    XmippIndexFlat2D(T* data_ptr, int column_count) : data(data_ptr), cols(column_count) {}

    T& operator[](int packed) const
    {
        return data[static_cast<std::ptrdiff_t>(packed) - cols - 1];
    }

    explicit operator bool() const { return data != NULL; }
};

template <class T>
class XmippIndexMatrix
{
public:
    T** rows;
    int row_lower;
    int col_lower;

    XmippIndexMatrix() : rows(NULL), row_lower(0), col_lower(0) {}

    class Row
    {
    public:
        T* data;
        int col_lower;
        T& operator[](int j) const
        {
            return data[static_cast<std::ptrdiff_t>(j) - col_lower];
        }
    };

    Row operator[](int i) const
    {
        Row row;
        row.data = rows[static_cast<std::ptrdiff_t>(i) - row_lower];
        row.col_lower = col_lower;
        return row;
    }

    explicit operator bool() const { return rows != NULL; }
    bool operator==(std::nullptr_t) const { return rows == NULL; }
    bool operator!=(std::nullptr_t) const { return rows != NULL; }
};

template <class T>
class XmippIndexVolume
{
public:
    T*** slices;
    int slice_lower;
    int row_lower;
    int col_lower;

    XmippIndexVolume() : slices(NULL), slice_lower(0), row_lower(0), col_lower(0) {}

    class Row
    {
    public:
        T* data;
        int col_lower;
        T& operator[](int j) const
        {
            return data[static_cast<std::ptrdiff_t>(j) - col_lower];
        }
    };

    class Slice
    {
    public:
        T** rows;
        int row_lower;
        int col_lower;
        Row operator[](int i) const
        {
            Row row;
            row.data = rows[static_cast<std::ptrdiff_t>(i) - row_lower];
            row.col_lower = col_lower;
            return row;
        }
    };

    Slice operator[](int k) const
    {
        Slice slice;
        slice.rows = slices[static_cast<std::ptrdiff_t>(k) - slice_lower];
        slice.row_lower = row_lower;
        slice.col_lower = col_lower;
        return slice;
    }

    explicit operator bool() const { return slices != NULL; }
    bool operator==(std::nullptr_t) const { return slices == NULL; }
    bool operator!=(std::nullptr_t) const { return slices != NULL; }
};

/* Memory managing --------------------------------------------------------- */
///@defgroup MemoryManaging Memory management for numerical recipes
/// @ingroup DataLibrary
//@{
/** Ask memory for any type vector.
    The valid values range from v[nl] to v[nh]. If no memory is available
    an exception is thrown. NULL is returned if nh is not greater than nl*/
template <class T> void ask_Tvector(XmippIndexPtr<T> &v, int nl, int nh)
{
    const std::size_t count = xmipp_index_count(nl, nh);
    if (count > 1)
    {
        T* base = static_cast<T*>(malloc(count * sizeof(T)));
        if (!base) REPORT_ERROR("allocation failure in vector()");
        v = XmippIndexPtr<T>(base, nl);
    }
    else v = XmippIndexPtr<T>();
}

/** Free memory associated to any type vector.
    After freeing v=NULL*/
template <class T> void free_Tvector(XmippIndexPtr<T> &v, int, int)
{
    if (v.data != NULL)
    {
        free(static_cast<void*>(v.data));
        v = XmippIndexPtr<T>();
    }
}

/** Ask memory for any type matrix.
    The valid values range from v[nrl][ncl] to v[nrh][nch].
    If no memory is available an exception is thrown. NULL is returned if any
    nh is not greater than its nl*/
template <class T> void ask_Tmatrix(XmippIndexMatrix<T> &m, int nrl, int nrh,
                                    int ncl, int nch)
{
    const std::size_t row_count = xmipp_index_count(nrl, nrh);
    const std::size_t column_count = xmipp_index_count(ncl, nch);
    if (row_count > 1 && column_count > 1)
    {
        T** row_base = static_cast<T**>(malloc(row_count * sizeof(T*)));
        if (!row_base) REPORT_ERROR( "allocation failure 1 in matrix()");
        for (int i = nrl;i <= nrh;i++)
        {
            T* base = static_cast<T*>(malloc(column_count * sizeof(T)));
            if (!base) REPORT_ERROR( "allocation failure 2 in matrix()");
            row_base[static_cast<std::ptrdiff_t>(i) - nrl] = base;
        }
        m.rows = row_base;
        m.row_lower = nrl;
        m.col_lower = ncl;
    }
    else m = XmippIndexMatrix<T>();
}

/** Free memory associated to any type matrix.
    After freeing v=NULL*/
template <class T> void free_Tmatrix(XmippIndexMatrix<T> &m, int nrl, int nrh,
                                     int, int)
{
    if (m.rows != NULL)
    {
        for (int i = nrh;i >= nrl;i--)
            free(static_cast<void*>(m.rows[static_cast<std::ptrdiff_t>(i) - nrl]));
        free(static_cast<void*>(m.rows));
        m = XmippIndexMatrix<T>();
    }
}

/** Ask memory for any type voliume.
    The valid values range from v[nsl][nrl][ncl] to v[nsh][nrh][nch].
    If no memory is available an exception is thrown. NULL is returned if any
    nh is not greater than its nl. */
template <class T> void ask_Tvolume(XmippIndexVolume<T> &m, int nsl, int nsh, int nrl,
                                    int nrh, int ncl, int nch)
{
    const std::size_t slice_count = xmipp_index_count(nsl, nsh);
    const std::size_t row_count = xmipp_index_count(nrl, nrh);
    const std::size_t column_count = xmipp_index_count(ncl, nch);
    if (slice_count > 1 && row_count > 1 && column_count > 1)
    {
        T*** slice_base = static_cast<T***>(malloc(slice_count * sizeof(T**)));
        if (!slice_base) REPORT_ERROR( "allocation failure 1 in matrix()");
        for (int k = nsl;k <= nsh;k++)
        {
            T** row_base = static_cast<T**>(malloc(row_count * sizeof(T*)));
            if (!row_base) REPORT_ERROR( "allocation failure 2 in matrix()");
            for (int i = nrl;i <= nrh;i++)
            {
                T* base = static_cast<T*>(malloc(column_count * sizeof(T)));
                if (!base) REPORT_ERROR( "allocation failure 2 in matrix()");
                row_base[static_cast<std::ptrdiff_t>(i) - nrl] = base;
            }
            slice_base[static_cast<std::ptrdiff_t>(k) - nsl] = row_base;
        }
        m.slices = slice_base;
        m.slice_lower = nsl;
        m.row_lower = nrl;
        m.col_lower = ncl;
    }
    else m = XmippIndexVolume<T>();
}

/** Free memory associated to any type volume.
    After freeing v=NULL*/
template <class T> void free_Tvolume(XmippIndexVolume<T> &m, int nsl, int nsh,
                                     int nrl, int nrh, int, int)
{
    if (m.slices != NULL)
    {
        for (int k = nsh;k >= nsl;k--)
        {
            T** rows = m.slices[static_cast<std::ptrdiff_t>(k) - nsl];
            for (int i = nrh;i >= nrl;i--)
                free(static_cast<void*>(rows[static_cast<std::ptrdiff_t>(i) - nrl]));
            free(static_cast<void*>(rows));
        }
        free(static_cast<void*>(m.slices));
        m = XmippIndexVolume<T>();
    }
}

/** Allocates memory.
 * Adapted from Bsofts bfree
 *
 * It is called exactly like malloc, with the following enhancements:
 *
 * - If allocation of zero bytes are requested it notifies the user.
 * - NO LONGER TRUE: Successfully allocated memory is zeroed
 * - Allocation is attempted and an error message is printed on failure.
 * - All failures return a NULL pointer to allow error handling from
 *    calling functions.
 *
 * returns char* : a pointer to the memory (NULL on failure)
 */
char* askMemory(unsigned long size);

/** Frees allocated memory.
 * Adapted from Bsofts bfree
 *
 * It is called exactly like free, with the following enhancements:
 *  - If freeing fails an error message is printed.
 *  - the pointer is reset to NULL
 *
 * returns int: 0 = success, -1 = failure.
*/
int freeMemory(void* ptr, unsigned long memsize);

//@}
#endif

