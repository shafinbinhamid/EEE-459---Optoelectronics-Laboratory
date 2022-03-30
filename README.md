# Automated-Hardware-Chess-Machine

This repository contains all codes for the 'Automated Hardware Chess Machine' project made for EEE-211 Numerical Techniques and Analysis Laboratory.

Inside @Board folder, you can find the logic designed to detect any illegal moves made by the human player.

Inside Control folder, there is control logic that is used to move the pieces around using an external control system including an arduino microcontroller with some brushless DC motors.

Inside Workdir, the codes are meant to synchronize everything. That is, first the game waits for the human player to make a move; then after a move has been made, the game checks if the move is legal; if legal, it sends the updated board picture to a chess machine called Stockfish from where the computer's move is generated and finally the control logic moves a particular piece according to the computer's command. Then, the game waits for the human player's next move. 

## Team

- Shafin Bin Hamid (github - https://github.com/shafinbinhamid) 
- Mir Sayeed Mohammad (github - https://github.com/ClockWorkKid)
 