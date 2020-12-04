include("../Network_qse_run1/Dependencies.jl")
gr()

a = Network_qse.extract_partition_function()



qse_alt = load("../Network_qse_run1/qse_tables/vary_R/QSE_table2.jld")["data"]
nse_alt = load("../Network_qse_run1/qse_tables/vary_R/NSE_table2.jld")["data"]

params_qse_alt = load("../Network_qse_run1/qse_tables/vary_R/QSE_params2.jld")
params_nse_alt = load("../Network_qse_run1/qse_tables/vary_R/NSE_params2.jld")

qse = load("../Network_qse_run1/NSE_higher_QSE/higher_4dec/QSE_table.jld")["data"]
nse = load("../Network_qse_run1/NSE_higher_QSE/higher/NSE_table.jld")["data"]

params_qse = load("../Network_qse_run1/NSE_higher_QSE/higher_4dec/QSE_params.jld")
params_nse = load("../Network_qse_run1/NSE_higher_QSE/higher/NSE_params.jld")


cl_qse = params_qse["srange"]
cl_nse = params_nse["srange"]
trange = params_qse["trange"]
yrange = params_qse["yrange"]
rrange = params_qse["rrange"]

AexNse(temp::Int64) = sum([nse[i,1,temp,1] * a[i].A for i in 1:size(a,1)])
AexQse(temp::Int64, r) = sum([qse[i,1,temp,1,r] * a[i].A for i in 1:size(a,1)])
ZexNse(temp::Int64) = sum([nse[i,1,temp,1] * a[i].Z for i in 1:size(a,1)])
ZexQse(temp::Int64, r) = sum([qse[i,1,temp,1,r] * a[i].Z for i in 1:size(a,1)])


p = Plots.plot(trange, cl_nse, label = "NSE", axis = "T₉", yaxis = "Clustersize", xaxis=:flip,legend = :bottomright)
hline!([cl_qse[1], cl_qse[end]], fillrange=[cl_qse[1] cl_qse[end]], alpha = 0.2,label = "QSE", axis = "T₉", yaxis = "Clustersize", xaxis=:flip)
Plots.savefig(p,"../Network_qse_run1/NSE_higher_QSE/higher_4dec/clustersize_new.pdf")


p = Plots.plot(trange, AexNse.(1:size(trange, 1)), label = "<A> NSE", linestyle = :dash, c = :blue,axis = "T₉", yaxis = "<Z>, <A>", xaxis=:flip, legend=:topleft)
plot!(trange, ZexNse.(1:size(trange, 1)), label = "<Z> NSE", linestyle = :dash, c = :red, axis = "T₉", yaxis = "<Z>, <A>", xaxis=:flip, legend=:topleft)
plot!(trange, AexQse.(1:size(trange, 1), size(cl_qse,1)), label = "<A> QSE", c = :blue, axis = "T₉", yaxis = "<Z>, <A>", xaxis=:flip, legend=:topleft)
plot!(trange, ZexQse.(1:size(trange, 1), size(cl_qse,1)), label = "<Z> QSE", c = :red, axis = "T₉", yaxis = "<Z>, <A>", xaxis=:flip, legend=:topleft)
plot!([trange trange], [ZexQse.(1:size(trange, 1), 1) ZexQse.(1:size(trange, 1), size(cl_qse,1))],  label = :false, fillrange=[ZexQse.(1:size(trange, 1), size(cl_qse,1)) ZexQse.(1:size(trange, 1), 1)], fillalpha=:0.1, color =:red)
plot!([trange trange], [AexQse.(1:size(trange, 1), 1) AexQse.(1:size(trange, 1), size(cl_qse,1))],  label = :false, fillrange=[AexQse.(1:size(trange, 1), size(cl_qse,1)) AexQse.(1:size(trange, 1), 1)], fillalpha=:0.1, color =:blue)

Plots.savefig("../Network_qse_run1/NSE_higher_QSE/higher_4dec/compare_below.png")



plot(trange, qse[Network_qse.find_el("Ni56", a), 1,:,1,1], xlabel = "T₉", ylabel = "Xᵢ", label = "QSE Ni62", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
plot!(trange, nse[Network_qse.find_el("Ni56", a), 1,:,1], xlabel = "T₉", ylabel = "Xᵢ", label = "NSE Ni62", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
savefig("Ni56.png")

plot!(trange, qse[Network_qse.find_el("H1", a), 1,:,1,1], xlabel = "T₉", ylabel = "Xᵢ", label = "QSE p", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
plot!(trange, nse[Network_qse.find_el("H1", a), 1,:,1], xlabel = "T₉", ylabel = "Xᵢ", label = "NSE p", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
savefig("H1.png")


plot!(trange, qse[Network_qse.find_el("1n", a), 1,:,1,1], xlabel = "T₉", ylabel = "Xᵢ", label = "QSE n", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
plot!(trange, nse[Network_qse.find_el("1n", a), 1,:,1], xlabel = "T₉", ylabel = "Xᵢ", label = "NSE n", xaxis =:flip, yaxis = :log, ylims = (1e-10, 1))
savefig("compare_all.png")




function species_str(a_i)
    if a_i.name == "1n"
        x = "n+p"
        y = " "
    elseif a_i.name == "H1"
        x = "n+p"
        y = " "
    elseif a_i.name == "H2"
        x = "H"
        y = "2"
    elseif a_i.name == "H3"
        x = "H"
        y = "3"
    elseif a_i.name == "He3"
        x = "He"
        y = "3"
    elseif a_i.name == "He4"
        x = "He"
        y = "4"
    elseif length(a_i.name) == 4.0
        y = a_i.name[3:4]
        x = a_i.name[1:2]
    else
        y = a_i.name[2:3]
        x = a_i.name[1]
    end
    return x, y
end


mutable struct NSE end

@recipe function f(::NSE, n::Integer, el::Integer, range::Vector, xaxis::String,
                    yVal::Array, customcolor = :green; add_marker = true)
    linecolor   --> customcolor
    seriestype  :=  :path
    markershape --> (add_marker ? :star5 : :none)
    #ytickfont -> font(20, "Courier")
    #markercolor --> "green"
    label --> false
    #xguide --> "Yₑ"
    xguide --> xaxis
    yguide --> "Xᵢ"
    #xtickfont --> font(50, "Courier")
    #ytickfont --> font(50, "Courier")
    yaxis --> :log
    linewidth = 4
    ylims --> (1e-2, 1)
    size --> (1920,1080)
    #linewidth = 4
    delete!(plotattributes, :true)
    #trange, all[el,n,:,1] # particle - Y - : - rho
    #yrange, all[el,:,n,1] # particle - : - tem - rho
    range, yVal
end



#col = shuffle(shuffle(Base.range(colorant"green", stop=colorant"red", length=size(a,1))))
#col = cgrad(:roma, rev = true, size(a,1), categorical = false, scale = :exp)
#col = cgrad(:tab20c, size(a,1), categorical = true, scale = :linear)
#col = shuffle(shuffle(Base.range(HSV(0,1,1), stop=HSV(-360,1,1),length=size(a,1))))
#col = shuffle(shuffle(Base.range(HSV(0,1,1), stop=HSV(-360,1,1),length=size(a,1))))

col = distinguishable_colors(size(a,1), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
shuffle!(col)


cl_qse

animY = @animate for i ∈ reverse(1:size(cl_qse,1)) # timeframes T
    println(i)
    Plots.plot()

    for (k, el) in zip(Iterators.countfrom(3), a[3:end])#for (k, el) in enumerate(a[3:end])
        #println(k, el)
        ind = argmax(filter!(!isnan, nse[k,1,:,1]))
        if nse[k,1,:,1][ind] > 0.01
            plot!(trange, nse[k, 1, :, 1], label = nothing, color = "grey", yaxis = :log)
            x,y = species_str(a[k])

            annotate!(trange[ind], 1.2*nse[k,1,:,1][ind], Plots.text(L"{}^{%$y}\!\textrm{%$x}", 25, "grey"))
        end
    end
    ind = argmax(filter!(!isnan, qse[1,1,:,1,i]))
    ind = argmax(filter!(!isnan, qse[1,1,:,1,i]))
    Plots.plot!(trange, qse[1,1,:,1, i]+qse[2,1,:,1, i], c = "red", label = "QSE p + n",ls=:dashdot,legendfontsize=12)
    Plots.plot!(trange, nse[1,1,:,1]+nse[2,1,:,1], c = "grey", label = "NSE p + n", ls=:dashdot,legendfontsize=12)
    annotate!(trange[50], 1.2, Plots.text(L"\mathrm{X}_{\mathrm{Cl}} \;\;\;\;\;\;\; = ", 19, "black"))
    annotate!(trange[50], 0.9, Plots.text(L"\overline{\mathrm{X}}_{\mathrm{NSE}} \;\;\;\;\; = ", 19, "black"))
    annotate!(trange[50], 0.7, Plots.text(L"\mathrm{X}_\mathrm{Cl}/\overline{\mathrm{X}}_\mathrm{NSE} = ", 19, "black"))
    annotate!(trange[40], 1.2, Plots.text(lpad(rpad(string(round(cl_qse[i]; sigdigits = 2, base = 10)),4), 4), 19, "black"))
    annotate!(trange[40], 0.9, Plots.text(lpad(rpad(string(round(mean(cl_nse[i]); sigdigits = 2, base = 10)),4), 4), 19, "black"))
    annotate!(trange[40], 0.7, Plots.text(lpad(rpad(string(round(cl_qse[i]/mean(cl_nse[i]); sigdigits = 2, base = 10)),4), 4), 19, "black"))

    annotate!(trange[90], 1.2, Plots.text(string(round(rrange[1]; sigdigits = 2, base = 10))*" g/cm³", 19, "black"))
    #annotate!(trange[end - 5], 1.2*(qse[1,1,:,1, i]+qse[2,1,:,1, i]), Plots.text(L"n + p", 15, "red"))
    #annotate!(trange[end - 5], 1.2*(nse[1,1,:,1]+nse[2,1,:,1]),Plots.text(L"n + p", 15, "grey"))
    for (k, el) in zip(Iterators.countfrom(3), a[3:end]) #enumerate(a[3:end])
        ind = argmax(filter!(!isnan, qse[k,1,:,1,i]))
        if qse[k,1,:,1,i][ind] > 0.01
            #plot!(nse,i, k, yrange, "Yₑ", all[k,:,i,1], colT[rand(1:30)], marker = (:circle,2),
            Plots.plot!(trange, qse[k,1,:,1, i],
            yaxis = :log,
            lw = 4,
            label = false,
            ylims = (1e-2, 1.5),
            yticks = ([1e-3, 1e-2, 1e-1, 1], ["10⁻³", "10⁻²", "10⁻¹", "1"]),
            size = (1920,1380),
            xlabel = "T₉ [K]",
            ylabel = "Xᵢ",
            c = col[k],
            #marker = (:circle,2),
            #markercolor = "white",
            linewidth = 2,
            guidefont=font(23),
            xtickfont = font(16),
            ytickfont = font(16),
            thickness_scaling = 1.3,
            margin=5Plots.mm)
            x,y = species_str(a[k])
            if qse[k,1,:,1, i][ind] > 0.001 && qse[k,1,:,1, i][ind] < 0.01
                annotate!(trange[ind], 1.2*qse[k,1,:,1, i][ind],
                Plots.text(L"{}^{%$y}\!\textrm{%$x}", 15, col[k]))
            else
                annotate!(trange[ind], 1.2*qse[k,1,:,1, i][ind],
                Plots.text(L"{}^{%$y}\!\textrm{%$x}", 20, col[k]))
            end
        end
    end
end

gif(animY, "x_cl_evol_2.mp4", fps = 12)
