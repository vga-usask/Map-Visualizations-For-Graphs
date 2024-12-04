#ifndef RPLayoutAlgorithm_hpp
#define RPLayoutAlgorithm_hpp

#include "RPGraphLayout.hpp"

namespace RPGraph
{
    class LayoutAlgorithm
    {
    public:
        LayoutAlgorithm(GraphLayout &layout);
        ~LayoutAlgorithm();
        GraphLayout &layout;

        virtual bool sync_layout() =0; // write current layout to `layout'.
        virtual void sync_layout_void() =0; // write current layout to `layout'.
    };
}

#endif

