
#pragma once

#include <algorithm>
#include <cmath>

namespace Blur {

    void std_to_box(int boxes[], float sigma, int n)  
    {
        // ideal filter width
        float wi = std::sqrt((12*sigma*sigma/n)+1); 
        int wl = std::floor(wi);  
        if(wl%2==0) wl--;
        int wu = wl+2;
                    
        float mi = (12*sigma*sigma - n*wl*wl - 4*n*wl - 3*n)/(-4*wl - 4);
        int m = std::round(mi);
                    
        for(int i=0; i<n; i++) 
            boxes[i] = ((i < m ? wl : wu) - 1) / 2;
    }

    //!
    //! \fn void horizontal_blur(float * in, float * out, int w, int h, int r)    
    //!
    //! \brief this function performs the horizontal blur pass for box blur. 
    //!
    //! \param[in,out] in       source channel
    //! \param[in,out] out      target channel
    //! \param[in] w            image width
    //! \param[in] h            image height
    //! \param[in] r            box dimension
    //!
    void horizontal_blur(float * in, float * out, int w, int h, int r) 
    {
        float iarr = 1.f / (r+r+1);
        #pragma omp parallel for
        for(int i=0; i<h; i++) 
        {
            int ti = i*w, li = ti, ri = ti+r;
            float fv = in[ti], lv = in[ti+w-1], val = (r+1)*fv;

            for(int j=0; j<r; j++) val += in[ti+j];
            for(int j=0  ; j<=r ; j++) { val += in[ri++] - fv      ; out[ti++] = val*iarr; }
            for(int j=r+1; j<w-r; j++) { val += in[ri++] - in[li++]; out[ti++] = val*iarr; }
            for(int j=w-r; j<w  ; j++) { val += lv       - in[li++]; out[ti++] = val*iarr; }
        }
    }

    //!
    //! \fn void total_blur(float * in, float * out, int w, int h, int r)   
    //!
    //! \brief this function performs the total blur pass for box blur. 
    //!
    //! \param[in,out] in       source channel
    //! \param[in,out] out      target channel
    //! \param[in] w            image width
    //! \param[in] h            image height
    //! \param[in] r            box dimension
    //!
    void total_blur(float * in, float * out, int w, int h, int r) 
    {
        float iarr = 1.f / (r+r+1);
        #pragma omp parallel for
        for(int i=0; i<w; i++) 
        {
            int ti = i, li = ti, ri = ti+r*w;
            float fv = in[ti], lv = in[ti+w*(h-1)], val = (r+1)*fv;
            for(int j=0; j<r; j++) val += in[ti+j*w];
            for(int j=0  ; j<=r ; j++) { val += in[ri] - fv    ; out[ti] = val*iarr; ri+=w; ti+=w; }
            for(int j=r+1; j<h-r; j++) { val += in[ri] - in[li]; out[ti] = val*iarr; li+=w; ri+=w; ti+=w; }
            for(int j=h-r; j<h  ; j++) { val += lv     - in[li]; out[ti] = val*iarr; li+=w; ti+=w; }
        }
    }

    //!
    //! \fn void box_blur(float * in, float * out, int w, int h, int r)    
    //!
    //! \brief this function performs a box blur pass. 
    //!
    //! \param[in,out] in       source channel
    //! \param[in,out] out      target channel
    //! \param[in] w            image width
    //! \param[in] h            image height
    //! \param[in] r            box dimension
    //!
    void box_blur(float *& in, float *& out, int w, int h, int r) 
    {
        std::swap(in, out);
        horizontal_blur(out, in, w, h, r);
        total_blur(in, out, w, h, r);
        // Note to myself : 
        // here we could go anisotropic with different radiis rx,ry in HBlur and TBlur
    }

    //!
    //! \fn void fast_gaussian_blur(float * in, float * out, int w, int h, float sigma)   
    //!
    //! \brief this function performs a fast Gaussian blur. Applying several
    //! times box blur tends towards a true Gaussian blur. Three passes are sufficient
    //! for good results. For further details please refer to :  
    //! http://blog.ivank.net/fastest-gaussian-blur.html
    //!
    //! \param[in,out] in       source channel
    //! \param[in,out] out      target channel
    //! \param[in] w            image width
    //! \param[in] h            image height
    //! \param[in] r            box dimension
    //!
    void fast_gaussian_blur(float *& in, float *& out, int w, int h, float sigma) 
    {
        // sigma conversion to box dimensions
        int boxes[3];
        std_to_box(boxes, sigma, 3);
        box_blur(in, out, w, h, boxes[0]);
        box_blur(out, in, w, h, boxes[1]);
        box_blur(in, out, w, h, boxes[2]);
    }
};