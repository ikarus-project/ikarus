//
// Created by Alex on 19.05.2021.
//

#ifndef IKARUS_NODEDATAINTERFACE_H
#define IKARUS_NODEDATAINTERFACE_H
#include <concepts>

class NodeDataInterface {

};

namespace Ikarus {
    template <typename NodeDataInterfaceType>
    concept NodeData = std::As requires (NodeDataInterfaceType node ){
        typename Node::ctype;
        {   node.addDof(grid) } ->  std::same_as<typename FEConceptType::MatrixType>;
    };

};


#endif //IKARUS_NODEDATAINTERFACE_H
