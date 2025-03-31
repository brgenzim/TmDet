// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>

/**
 * @brief namespace for tmdet system
 *
 * @namespace Tmdet
 * @namespace System
 */
namespace Tmdet::System {

    /**
     * @brief simple class for console progress bar
     */
    class ProgressBar {
        private:
            /**
             * @brief  
             */
            std::string _title;

            /**
             * @brief size of the progress bar in character
             */
            unsigned int _width = 80;

            /**
             * @brief number of ticks
             */
            unsigned long _numTicks;

            /**
             * @brief current tick
             */
            unsigned long _currentTick = 0;

            /**
             * @brief last tick when bar was shown
             */
            unsigned long _lastTick = 0;

            /**
             * @brief progress
             */
            double _progress;
            
            /**
             * @brief doneCharacter
             */
            char _doneCharacter = '#';

            /**
             * @brief todoCharacter
             */
            char _todoCharacter = '.';

            /**
             * @brief display percentage flag
             */
            bool _displayPercentage = false;

            /**
             * @brief display tasks done flag
             */
            bool _displayTaskDone = false;

            /**
             * @brief display the progress bar
             */
            void showBar();
            

        public:
            
            /**
             * @brief Set the Title of the progress bar
             * 
             * @param title 
             */
            void setTitle(const char* title) {
                _title = title;
            }

            /**
             * @brief Set the width of the progress bar (in character)
             * 
             * @param width 
             */
            void setWidth(const unsigned int width) {
                _width = width;
            }

            /**
             * @brief Set the number of ticks
             * 
             * @param numTicks 
             */
            void setNumTicks(const unsigned long numTicks) {
                _numTicks = numTicks;
            }

            /**
             * @brief Set the Done Character
             * 
             * @param doneCharacter 
             */
            void setDoneCharacter(const char doneCharacter) {
                _doneCharacter = doneCharacter;
            }

            /**
             * @brief Set the Todo Character
             * 
             * @param todoCharacter 
             */
            void setTodoCharacter(const char todoCharacter) {
                _todoCharacter = todoCharacter;
            }
            
            /**
             * @brief increase tick by one
             */
            void tick();
            
            /**
             * @brief set display percentage flag
             */
            void displayPercentage() {
                _displayPercentage = true;
            }
        
            /**
             * @brief set display task done flag
             */
            void displayTasksDone() {
                _displayTaskDone = true;
            }

            /**
             * @brief end of progress bar
             */
            void end() {
                std::cout << "\x1b[2K";
                _currentTick = _numTicks;
                _displayPercentage = false;
                _displayTaskDone = false;
                showBar();
                std::cout << " Done" << std::endl;
            }
            
    };
}
