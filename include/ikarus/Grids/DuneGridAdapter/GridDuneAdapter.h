//
// Created by Alex on 25.05.2021.
//

#pragma once

//Adpator Pattern 4.1 GOF Book
template<typename DuneGridType> //TODO: define Concept
class GridDuneAdapter {
    template<class... Args>
    GridDuneAdapter(Args... args)
    :underLyingDuneGrid(args...)
    {}
private:
    DuneGridType underLyingDuneGrid;
};