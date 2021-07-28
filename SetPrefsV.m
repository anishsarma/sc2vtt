        switch allvcases{vInd}
            case 'Ideal'
                popStructBool = 0;
                ifnDelayedBool = 0;
                skipRate = 0;
            case 'Host-Varying'
                popStructBool = 0;
                ifnDelayedBool = 0;
                skipRate = 1;
            case 'Time-Varying'
                popStructBool = 0;
                ifnDelayedBool = 1;
                skipRate = 0;
            case 'Time- and Host-Varying'
                popStructBool = 0;
                ifnDelayedBool = 1;
                skipRate = 1;
            case 'Time- and Host-Varying with Population Stucture'
                popStructBool = 1;
                ifnDelayedBool = 1;
                skipRate = 1;
        end
