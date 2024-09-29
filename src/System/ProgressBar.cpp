#include <System/ProgressBar.hpp>

namespace Tmdet::System {
    
    void ProgressBar::tick() {
        _currentTick++;
        if ((double)(_currentTick - _lastTick) * _width / _numTicks > 0.8 ) {
            showBar();
            _lastTick = _currentTick;
        }
    }

    void ProgressBar::showBar() {
        _progress = (double)_currentTick / _numTicks;
        std::cout << '\r' << _title << "[";
        for(unsigned int i = 0; i < _width; i++) {
            std::cout << (i<_progress*_width?_doneCharacter:_todoCharacter);
        }
        std::cout << "] ";
        if (_displayPercentage) {
            std::cout << (int)(100*_progress) << "%";
        }
        if (_displayTaskDone ) {
            std::cout << " (" << _currentTick << '/' << _numTicks << ')';
        }
        std::cout << std::flush;
    }

}