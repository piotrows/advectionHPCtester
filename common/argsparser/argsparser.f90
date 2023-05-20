MODULE ArgsParser
    IMPLICIT NONE
CONTAINS
    
    LOGICAL FUNCTION CheckForHelp()
        CHARACTER(len=128)  :: arg
        INTEGER             :: i
         
        CheckForHelp = .FALSE.
        DO i = 1, command_argument_count() 
            CALL getarg(i, arg)
            IF (TRIM(arg) == "-h") THEN
                CheckForHelp = .TRUE.
            END IF
        END DO
        RETURN
    END FUNCTION CheckForHelp

    SUBROUTINE ParseArgument(ArgName, ArgType, ArgValue, SuccessFlag, StatusCorrect, StatusIncorrect, UseDefault)
        CHARACTER(len=*)    :: ArgName
        INTEGER             :: ArgType
        CHARACTER(len=128)  :: ArgValue   
        INTEGER             :: SuccessFlag
        INTEGER             :: StatusCorrect
        INTEGER             :: StatusIncorrect
        LOGICAL             :: UseDefault

        INTEGER             :: i
        INTEGER             :: ArgCount
        
        CHARACTER(len=128)  :: arg
        


        ArgCount = command_argument_count() 
        SuccessFlag = 0
        IF (ArgType /= 10) THEN
            DO i = 1, ArgCount
                CALL getarg(i, arg)
                ! WRITE (*,*) arg
                IF (TRIM(ArgName) == TRIM(arg)) THEN
                    IF (i+1 <= ArgCount) THEN
                        CALL getarg(i+1, ArgValue)
                        SuccessFlag = SuccessFlag + 1
                    END IF
                END IF
            END DO
        ELSE
            DO i = 1, ArgCount
                CALL getarg(i, arg)
                ! WRITE (*,*) arg
                IF (TRIM(ArgName) == TRIM(arg)) THEN
                    SuccessFlag = SuccessFlag + 1
                    ArgValue = "1"//REPEAT(" ",127)
                END IF
            END DO
        END IF

        IF (SuccessFlag == 1) THEN
            StatusCorrect = StatusCorrect + 1
        ELSE
            ArgValue = ""
            IF (SuccessFlag == 0) THEN
                IF (UseDefault .eqv. .FALSE.) THEN
                    StatusIncorrect = StatusIncorrect + 1
                    write(*,*) "ERROR: ", TRIM(ArgName), " not found in argument list!"
                END IF
                
            END IF

            IF (SuccessFlag > 1) THEN
                write(*,*) "ERROR: ", TRIM(ArgName), " used multiple times!"
                StatusIncorrect = StatusIncorrect + 1
            END IF
        END IF
    END SUBROUTINE ParseArgument

    INTEGER FUNCTION ParseArgumentInt(ArgName, UseDefault, Default, StatusCorrect, StatusIncorrect, PrintHelp, Help)
        CHARACTER(len=*)    ArgName
        LOGICAL             UseDefault
        INTEGER             Default
        INTEGER             StatusCorrect
        INTEGER             StatusIncorrect
        LOGICAL             PrintHelp
        CHARACTER(len=*)    Help


        CHARACTER(len=128)  :: ArgValue
        integer             :: int
        integer             :: stat
        INTEGER             :: SuccessFlag

        IF (PrintHelp .eqv. .FALSE.) THEN   
            CALL ParseArgument(ArgName, 1, ArgValue, SuccessFlag, StatusCorrect, StatusIncorrect, UseDefault)
            IF (SuccessFlag == 1) THEN
                read(ArgValue,*,iostat=stat) int
                ParseArgumentInt = int
            ELSE
                IF (UseDefault .eqv. .TRUE.) THEN
                    ParseArgumentInt = Default

                ELSE
                    ParseArgumentInt = 0
                END IF
            END IF
        ELSE
            write(*,*) Help
        END IF

        RETURN
    END FUNCTION


    REAL FUNCTION ParseArgumentReal(ArgName, UseDefault, Default, StatusCorrect, StatusIncorrect, PrintHelp, Help)
        CHARACTER(len=*)    ArgName
        LOGICAL             UseDefault
        REAL                Default
        INTEGER             StatusCorrect
        INTEGER             StatusIncorrect
        LOGICAL             PrintHelp
        CHARACTER(len=*)    Help


        CHARACTER(len=128)  :: ArgValue
        REAL                :: int
        integer             :: stat
        INTEGER             :: SuccessFlag

        IF (PrintHelp .eqv. .FALSE.) THEN   
            CALL ParseArgument(ArgName, 1, ArgValue, SuccessFlag, StatusCorrect, StatusIncorrect, UseDefault)
            IF (SuccessFlag == 1) THEN
                read(ArgValue,*,iostat=stat) int
                ParseArgumentReal = int
            ELSE
                IF (UseDefault .eqv. .TRUE.) THEN
                    ParseArgumentReal = Default

                ELSE
                    ParseArgumentReal = 0
                END IF
            END IF
        ELSE
            write(*,*) Help
        END IF

        RETURN
    END FUNCTION



    LOGICAL FUNCTION ParseArgumentLogical(ArgName, UseDefault, Default, StatusCorrect, StatusIncorrect, PrintHelp, Help)
        CHARACTER(len=*)    ArgName
        LOGICAL             UseDefault
        LOGICAL             Default
        INTEGER             StatusCorrect
        INTEGER             StatusIncorrect
        LOGICAL             PrintHelp
        CHARACTER(len=*)    Help


        CHARACTER(len=128)  :: ArgValue
        INTEGER             :: SuccessFlag

        IF (PrintHelp .eqv. .FALSE.) THEN   
            CALL ParseArgument(ArgName, 10, ArgValue, SuccessFlag, StatusCorrect, StatusIncorrect, UseDefault)
            IF (SuccessFlag == 1) THEN
                IF (TRIM(ArgValue) == "1") THEN
                    ParseArgumentLogical = .TRUE.
                ELSE
                    ParseArgumentLogical = .FALSE.
                END IF
            ELSE
                IF (UseDefault .eqv. .TRUE.) THEN
                    ParseArgumentLogical = Default

                ELSE
                    ParseArgumentLogical = .FALSE.
                END IF
            END IF
        ELSE
            write(*,*) Help
        END IF

        RETURN
    END FUNCTION


    CHARACTER(len=128) FUNCTION ParseArgumentString(ArgName, UseDefault, Default, StatusCorrect, StatusIncorrect, PrintHelp, Help)
        CHARACTER(len=*)    ArgName
        LOGICAL             UseDefault
        CHARACTER(len=*)    Default
        INTEGER             StatusCorrect
        INTEGER             StatusIncorrect
        LOGICAL             PrintHelp
        CHARACTER(len=*)    Help


        CHARACTER(len=128)  :: ArgValue
        INTEGER             :: SuccessFlag

        IF (PrintHelp .eqv. .FALSE.) THEN   
            CALL ParseArgument(ArgName, 1, ArgValue, SuccessFlag, StatusCorrect, StatusIncorrect, UseDefault)
            IF (SuccessFlag == 1) THEN
                ! read(ArgValue,*,iostat=stat) int
                ParseArgumentString = ArgValue
            ELSE
                IF (UseDefault .eqv. .TRUE.) THEN
                    ParseArgumentString = Default

                ELSE
                    ParseArgumentString = repeat("",128)
                END IF
            END IF
        ELSE
            write(*,*) Help
        END IF

        RETURN
    END FUNCTION


END MODULE ArgsParser
