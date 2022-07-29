//
// Created by lex on 7/29/22.
//

#pragma once

namespace Ikarus {

#include <vector>
#include <ranges>

  /** \brief A multi-index
   */
  class MultiIndex
      : public std::vector<unsigned int>
  {
    unsigned int rangeOfEachComponent_;

  public:
    /** \brief Constructor with a given range for each digit
     * \param n Number of digits
     */
    MultiIndex(unsigned int numberOfComponents, unsigned int rangeOfEachComponent)
        : std::vector<unsigned int>(numberOfComponents),
          rangeOfEachComponent_(rangeOfEachComponent)
    {
      std::ranges::fill((*this), 0);
    }

    /** \brief Increment the MultiIndex */
    MultiIndex& operator++() {

      for (size_t i=0; i<size(); i++) {

        // Augment digit
        (*this)[i]++;

        // If there is no carry-over we can stop here
        if ((*this)[i]<rangeOfEachComponent_)
          break;

        (*this)[i] = 0;

      }
      return *this;
    }

    /** \brief Compute how many times you can call operator++ before getting to (0,...,0) again */
    size_t cycles() const {
      size_t result = 1;
      for (size_t i=0; i<this->size(); i++)
        result *= rangeOfEachComponent_;
      return result;
    }

  };


}  // namespace Ikarus


