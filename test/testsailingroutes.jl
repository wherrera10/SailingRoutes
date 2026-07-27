using Test
using SailingRoutes

function testweatherrouting()
    # Create polar data
    function create_sample_polar()
        winds = Float64[5, 10, 15, 20, 25]
        degrees = Float64[0, 30, 45, 60, 90, 120, 135, 150, 180]
        speeds = Float64[
            0.0  2.0  4.0  5.0  6.0;
            1.0  3.0  5.0  7.0  8.0;
            2.0  4.0  6.0  8.0  9.0;
            3.0  5.0  7.0  9.0  10.0;
            4.0  6.0  8.0  10.0 11.0;
            5.0  7.0  9.0  11.0 12.0;
            6.0  8.0  10.0 12.0 13.0;
            7.0  9.0  11.0 13.0 14.0;
            8.0  10.0 12.0 14.0 15.0
        ]
        return SailingPolar(winds, degrees, speeds)
    end

    # Define obstacles
    #=
        The data is selected so the best time path is slightly longer than the
        shortest length path. The forbidden regions are x, representing land or reef.
        The allowed sailing points are . and start and finish are S and F.
        Looking at this with the grid set as row 1 at the bottom, it is:

        x  .  .  F  .  .  x  .  x
        .  .  .  .  .  .  .  x  x
        x  .  .  x  x  x  .  .  .
        .  .  x  x  x  x  .  x  x
        x  .  .  .  x  x  .  x  .
        x  .  .  .  x  x  .  x  .
        .  .  .  .  x  .  .  x  .
        x  .  .  .  .  .  .  x  .
        .  .  .  S  .  x  .  .  .
    =#
    forbidden = falses(9, 9)  # Assuming a 9x9 grid for the example
    for idx in [[1, 6], [2, 1], [2, 8], [3, 5], [3, 8], [4, 1], [4, 5], [4, 6], [4, 8], [5, 1],
        [5, 5], [5, 6], [5, 8], [6, 3], [6, 4], [6, 5], [6, 6], [6, 8], [6, 9], [7, 1],
        [7, 4], [7, 5], [7, 6], [8, 8], [8, 9], [9, 1], [9, 7], [9, 9]]
        forbidden[idx[1], idx[2]] = true
    end

    # Create regional wind patterns
    function surfacebylongitude(lon)
        return lon < -155.04 ? SurfaceParameters(5, 2, 150, 0.5) :
               lon < -155.00 ? SurfaceParameters(45, 7, 150, 0.4) :
               SurfaceParameters(180.0, 10, 150, 0.3)
    end

    # Vary wind speeds over time
    function mutatetimeslices!(slices)
        for (i, slice) in enumerate(slices)
            for I in eachindex(slice)
                x = slice[I]
                slice[I] = GridPoint(
                    x.pt,
                    SurfaceParameters(x.sp.winddeg, x.sp.windkts * (1 + 0.01 * i),
                        x.sp.currentdeg, x.sp.currentkts),
                )
            end
        end
    end

    # Set up problem data
    startpos = SailingRoutes.GridPosition(1, 4)
    endpos = SailingRoutes.GridPosition(9, 4)
    pmat = [Position(19.78 - 1/60 + i/60, -155.0 - 5/60 + j/60) for i in 0:8, j in 0:8]
    gpoints = [SailingRoutes.GridPoint(pt, surfacebylongitude(pt.lon)) for pt in pmat]
    slices = [deepcopy(gpoints) for _ in 1:200]
    mutatetimeslices!(slices)

    routeprob = RoutingProblem(10.0, slices, forbidden, startpos, endpos)

    # Try to load polar data or use sample
    local sp, madepolar
    try
        sp = getpolardata("polar.csv")
        madepolar = false
    catch
        println("Using sample polar data (polar.csv not found)")
        sp = create_sample_polar()
        madepolar = true
    end

    println("Solving routing problem...")
    @time tp = minimumtimeroute(routeprob, sp; verbose = true)

    if isfinite(tp.duration)
        println("\nThe route taking the least time found was:")
        println("    Path: $(tp.path)")
        println("    Duration: $(floor(Int,tp.duration ÷ 60)) hours, $(round(Int, tp.duration % 60)) minutes")
        println("    Total points: $(length(tp.path))")
    else
        println("\nNo valid route found!")
    end

    @testset "Core utilities" begin

        @testset "Angle normalization" begin
            @test normalizeangledegrees(0.0) == 0.0
            @test normalizeangledegrees(360.0) == 0.0
            @test normalizeangledegrees(720.0) == 0.0
            @test normalizeangledegrees(-90.0) == 270.0
            @test normalizeangledegrees(450.0) == 90.0
            @test normalizeangledegrees(-123.4567) == 236.5433
        end

        @testset "angledifference" begin
            @test angledifference(350.0, 10.0) == 20.0
            @test angledifference(10.0, 350.0) == 20.0
            @test angledifference(0.0, 180.0) == 180.0
            @test angledifference(0.0, 0.0) == 0.0
            @test angledifference(90.0, 270.0) == 180.0
        end

        @testset "vmg" begin
            @test vmg(10.0, 0.0, 0.0) == 10.0
            @test isapprox(vmg(10.0, 0.0, 90.0), 0.0; atol = 1e-9)
            @test isapprox(vmg(10.0, 0.0, 180.0), -10.0; atol = 1e-9)
        end

        @testset "cartesian2polar / polar2cartesian round trip" begin
            for deg in (0.0, 45.0, 90.0, 200.0, 359.0)
                x, y = polar2cartesian(5.0, deg)
                r2, deg2 = cartesian2polar(x, y)
                @test isapprox(r2, 5.0; atol = 1e-9)
                @test isapprox(deg2, deg; atol = 1e-6)
            end
        end

        @testset "haversine" begin
            d0, _ = haversine(10.0, 20.0, 10.0, 20.0)
            @test isapprox(d0, 0.0; atol = 1e-6)

            _, bearing_north = haversine(0.0, 0.0, 1.0, 0.0)
            @test isapprox(bearing_north, 0.0; atol = 1e-6)

            _, bearing_east = haversine(0.0, 0.0, 0.0, 1.0)
            @test isapprox(bearing_east, 90.0; atol = 1e-6)
        end

        @testset "inverse_haversine round trip" begin
            lat2, lon2 = inverse_haversine(10.0, 20.0, 50_000.0, 45.0)
            d, bearing = haversine(10.0, 20.0, lat2, lon2)
            @test isapprox(d, 50_000.0; rtol = 1e-6)
            @test isapprox(bearing, 45.0; atol = 1e-6)
        end
        @testset "Longitude normalization" begin
            lat, lon = inverse_haversine(
                0.0,
                179.9,
                50_000.0,
                90.0,
            )
            @test -180 <= lon < 180
        end
        @testset "Bearing symmetry" begin
            d1, b1 = haversine(20, 30, 21, 31)
            d2, b2 = haversine(21, 31, 20, 30)
            @test isapprox(d1, d2; rtol = 1e-12)
            reversebearing = mod(b1 + 180, 360)
            @test isapprox(reversebearing, b2; atol = 2)
        end

        if !madepolar
            ep = 0.0005
            @testset "Speed and Angle Tests" begin
                @test abs((haversine(Position(0, 0), Position(0, 1))[1] - 111200) / 111200) < ep

                println("\n=== Testing Individual Functions ===")
                test_angle = 45.0
                test_speed = 15.0
                boat_speed = boatspeed(sp, test_angle, test_speed)
                @test abs(boat_speed - 0.07) < ep

                interp_speed = interpolatepolar(sp, test_angle, test_speed)
                @test abs(interp_speed - 0.07) < ep

                # Test best tack: 100 degrees with these parameters
                bestangle, bestvmg = findbesttack(sp, 45.0, test_speed)
                @show bestangle, bestvmg
                @test abs(bestangle - 100.0) < ep
                @test abs(bestvmg - 8.31685) < ep
            end
            @testset "Zero Wind Tests" begin
                zero_wind_sp = SailingPolar([0.0, 10.0], [0.0, 30.0], zeros(2, 2))
                test_angle = 45.0
                test_speed = 15.0
                boat_speed = boatspeed(zero_wind_sp, test_angle, test_speed)
                @test isapprox(boat_speed, 0.0)
                interp_speed = interpolatepolar(zero_wind_sp, test_angle, test_speed)
                @test isapprox(interp_speed, 0.0)
            end
        end

        @testset "Edge Cases" begin
            # Start == Finish, one-cell map
            begin
                pts = [GridPoint(Position(0.0, 0.0), SurfaceParameters());;]
                rp = RoutingProblem(
                    1.0,
                    [pts],
                    falses(1, 1),
                    (1, 1),
                    (1, 1),
                )
                tp = minimumtimeroute(rp, create_sample_polar())
                @test tp.duration == 0.0
                @test tp.path == [(1, 1)]
            end

            # Entire map blocked except start/finish
            begin
                pts = [GridPoint(Position(i, j), SurfaceParameters())
                       for i ∈ 1:3, j ∈ 1:3]
                obstacles = trues(3, 3)
                obstacles[1, 1] = false
                obstacles[3, 3] = false
                rp = RoutingProblem(
                    1.0,
                    TimeSlice[pts],
                    obstacles,
                    (1, 1),
                    (3, 3),
                )
                tp = minimumtimeroute(rp, create_sample_polar())
                @test !isfinite(tp.duration)
            end

            # Empty SailingPolar
            @test_throws ArgumentError SailingPolar(
                Float64[],
                Float64[],
                zeros(0, 0),
            )

            # Dateline crossing
            begin
                d, bearing = haversine(0.0, 179.9, 0.0, -179.9)
                @test d < 25_000        # should only be about 22 km
                @test isapprox(bearing, 90.0; atol = 0.5)
            end

            # North pole
            begin
                lat2, lon2 = inverse_haversine(
                    89.9,
                    0.0,
                    5000.0,
                    0.0,
                )
                @test lat2 <= 90.0
                @test lat2 > 89.9
            end

            # South pole
            begin
                lat2, lon2 = inverse_haversine(
                    -89.9,
                    0.0,
                    5000.0,
                    180.0,
                )
                @test lat2 >= -90.0
                @test lat2 < -89.9
            end

            # Zero current
            begin
                heading, sog = bestvectorspeed(
                    create_sample_polar(),
                    90.0,
                    180.0,
                    15.0,
                    0.0,
                    0.0,
                )
                @test sog > 0.0
                @test 0.0 <= heading <= 360.0
            end

            # Huge favorable current
            begin
                heading, sog = bestvectorspeed(
                    create_sample_polar(),
                    90.0,
                    0.0,
                    5.0,
                    90.0,
                    20.0,
                )
                @test sog > 20.0
            end

            # Huge adverse current
            begin
                heading, sog = bestvectorspeed(
                    create_sample_polar(),
                    90.0,
                    0.0,
                    5.0,
                    270.0,
                    20.0,
                )
                @test sog >= 0.0
            end
        end

        @testset "SailingPolar validation" begin
            @test_throws ArgumentError SailingPolar([5.0, 10.0], [0.0, 30.0], zeros(1, 2))
            @test_throws ArgumentError SailingPolar([-5.0, 10.0], [0.0, 30.0], zeros(2, 2))
            @test_throws ArgumentError SailingPolar([10.0, 5.0], [0.0, 30.0], zeros(2, 2))
        end

        @testset "TimedPath validation" begin
            @test_throws ArgumentError TimedPath(-1.0, [(1, 1)])
            @test_throws ArgumentError TimedPath(1.0, GridPosition[])
            tp_valid = TimedPath(5.0, [(1, 1), (2, 2)])
            @test tp_valid.duration == 5.0
            @test length(tp_valid.path) == 2
        end

        @testset "RoutingProblem validation" begin
            pts = [GridPoint(Position(0.0, 0.0), SurfaceParameters()) for i in 1:2, j in 1:2]
            timeframe = TimeSlice[pts]
            obstacles = falses(2, 2)

            @test_throws ArgumentError RoutingProblem(-1.0, timeframe, obstacles, (1, 1), (2, 2))
            @test_throws ArgumentError RoutingProblem(1.0, TimeSlice[], obstacles, (1, 1), (2, 2))
            @test_throws ArgumentError RoutingProblem(1.0, timeframe, obstacles, (1, 1), (5, 5))

            obstacles_blocked = falses(2, 2)
            obstacles_blocked[1, 1] = true
            @test_throws ArgumentError RoutingProblem(1.0, timeframe, obstacles_blocked, (1, 1), (2, 2))

            rp_valid = RoutingProblem(1.0, timeframe, obstacles, (1, 1), (2, 2))
            @test rp_valid.start == (1, 1)
            ok, _ = validateroutingproblem(rp_valid)
            @test ok
        end

        @testset "closestpoint" begin
            mat = [GridPoint(Position(Float64(i), Float64(j)), SurfaceParameters()) for i in 0:2, j in 0:2]
            cp = closestpoint(Position(1.1, 0.9), mat)
            @test cp.pt.lat == 1.0
            @test cp.pt.lon == 1.0
        end

        @testset "Time Slice Selection" begin
            rp = RoutingProblem(
                10.0,
                slices,
                forbidden,
                startpos,
                endpos,
            )
            @test gettimeslice(rp, 0.0) === rp.timeframe[1]
            @test gettimeslice(rp, 9.9) === rp.timeframe[1]
            @test gettimeslice(rp, 10.0) === rp.timeframe[2]
            @test gettimeslice(rp, 999999.0) === rp.timeframe[end]
        end

        @testset "route path sanity" begin
            if isfinite(tp.duration)
                @test tp.path[1] == startpos
                @test tp.path[end] == endpos
                @test length(tp.path) >= 2
            end
        end
    end

end

testweatherrouting()
