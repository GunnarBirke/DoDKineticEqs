function check_parameters(setup)
    eq_type = setup["eq_type"]
    fluxtype = setup["fluxtype"]
    J1_type = setup["J1_type"]
    J1_type_A = setup["J1_type_A"]

    tel_wave_fluxes = ["full", "altlr", "altrl", "central"]
    tel_wave_fluxes_todo = []
    transport_fluxes = ["upwind", "altlr", "altrl", "central"] # altlr und altrl könnten auch anstatt upwind geschrieben werden

    tel_wave_J1_types_altlr = ["upwind_classic", "upwind_symm", "upwind_diss_symm"]
    tel_wave_J1_types_central = ["central_classic", "central_symm"]
    tel_wave_J1_types_altrl = ["upwind_classic", "upwind_symm", "upwind_diss_symm"]
    tel_wave_J1_types_altlr_todo = []
    tel_wave_J1_types_central_todo = []
    tel_wave_J1_types_altrl_todo = []

    transport_J1_types_upwind = ["upwind_classic", "upwind_symm", "upwind_diss_symm"] #TODO die letzteren beiden müssen aber noch für die andere Richung implementiert werden
    transport_J1_types_central = ["central_classic", "central_symm"]
    transport_J1_types_upwind_todo = []
    transport_J1_types_central_todo = ["central_classic", "central_symm"]

    J1_type_A_fluxes = ["upwind_classic", "upwind_symm", "upwind_diss_symm"]
    J1_type_A_fluxes_todo = []

    ##################
    #Transport checks#
    ##################
    if eq_type == "transport"
        # alles außer full sollte hier funktionieren
        if !(fluxtype in transport_fluxes)
            throw(ArgumentError("Flusstyp ``$(fluxtype)`` not available for transport"))
        end
        if fluxtype == "altlr" || fluxtype == "altrl"
            setup["fluxtype"] = "upwind"
            fluxtype = "upwind"
            println("Fluxtype fixed to upwind, because alternating not possible for transport")
        end
        # J1-type checks für transport
        if fluxtype == "upwind" &&  !(J1_type in transport_J1_types_upwind) && J1_type in transport_J1_types_upwind_todo
            throw(ArgumentError("J1-Term noch nicht für upwind implementiert"))
        elseif fluxtype == "upwind" &&  !(J1_type in transport_J1_types_upwind) && !(J1_type in transport_J1_types_upwind_todo)
            throw(ArgumentError("J1-Term nicht für upwind möglich"))
        end
        if fluxtype == "central" &&  !(J1_type in transport_J1_types_central) && J1_type in transport_J1_types_central_todo
            throw(ArgumentError("J1-Term noch nicht für central implementiert"))
        elseif fluxtype == "central" &&  !(J1_type in transport_J1_types_upwind) && !(J1_type in transport_J1_types_central_todo)
            throw(ArgumentError("J1-Term nicht für central möglich"))
        end
    ###############################################
    #Telegraph/wave checks and therefore also heat#
    ###############################################
    elseif eq_type == "telegraph" || eq_type == "wave" || eq_type == "heat" # alle Fälle sollten beide Flüsse können
        if eq_type == "wave"
            println("Wahl von J1_type_A für Wellengleichung irrelevant: Term nicht vorhanden") 
        end
        if eq_type == "heat" && fluxtype == "full"
            throw(ArgumentError("Flusstyp ``full`` nicht für Wärmeleitungsgleichung implementiert!"))
        end
        if !(fluxtype in tel_wave_fluxes)
            if fluxtype in tel_wave_fluxes_todo
                throw(ArgumentError("Flusstyp ``$(fluxtype)`` nicht für Wellen bzw Telegraphengleichung implementiert!"))
            else
                throw(ArgumentError("Flusstyp ``$(fluxtype)`` nicht für Wellen bzw Telegraphengleichung möglich!"))
            end
        else
            #check J1-Konnektivität
            if fluxtype == "full" 
                println("Wahl von J1_type nicht relevant: Fluss eindeutig durch ``full`` festgelegt")
            end
            if fluxtype == "altlr" &&  !(J1_type in tel_wave_J1_types_altlr) && J1_type in tel_wave_J1_types_altlr_todo
                throw(ArgumentError("J1-Term noch nicht für altlr implementiert"))
            elseif fluxtype == "altlr" &&  !(J1_type in tel_wave_J1_types_altlr) && !(J1_type in tel_wave_J1_types_altlr_todo)
                throw(ArgumentError("J1-Term nicht für altlr möglich"))
            end
            if fluxtype == "altrl" &&  !(J1_type in tel_wave_J1_types_altrl) && J1_type in tel_wave_J1_types_altrl_todo
                throw(ArgumentError("J1-Term noch nicht für altrl implementiert"))
            elseif fluxtype == "altrl" &&  !(J1_type in tel_wave_J1_types_altrl) && !(J1_type in tel_wave_J1_types_altrl_todo)
                throw(ArgumentError("J1-Term nicht für altlr möglich"))
            end
            if fluxtype == "central" &&  !(J1_type in tel_wave_J1_types_central) && J1_type in tel_wave_J1_types_central_todo
                throw(ArgumentError("J1-Term noch nicht für central implementiert"))
            elseif fluxtype == "central" &&  !(J1_type in tel_wave_J1_types_central) && !(J1_type in tel_wave_J1_types_central_todo)
                throw(ArgumentError("J1-Term nicht für central möglich"))
            end
        end
    end
    ##################
    #J1_type_A checks#
    ##################
    if !(J1_type_A in J1_type_A_fluxes)
        if J1_type_A in J1_type_A_fluxes_todo
            throw(ArgumentError("J1_type_A Fluss ``$(J1_type_A)`` noch nicht implementiert"))
        else
            throw(ArgumentError("J1_type_A Fluss ``$(J1_type_A)`` nicht möglich!"))
        end
    end
    return setup
    # zu checken:
    # passt FLuss zu Gleichung?
    # passt J1_type zu Fluss?
    # zu J1_type_A eigentlich nichts zu checken
    # Exisiert Fluss zu Asymp-Limit Fluss? (z.B. keine Idee für full und aktuell nichts für central implementiert - sollte aber eigenltich aus Submatrix-Implementierung folgen)
end
