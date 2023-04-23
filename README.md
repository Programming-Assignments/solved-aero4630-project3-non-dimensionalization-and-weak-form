Download Link: https://assignmentchef.com/product/solved-aero4630-project3-non-dimensionalization-and-weak-form
<br>
<h1>Part 1a: Non-dimensionalization and Weak Form</h1>

The governing equation

div<em>σ </em>+ <em>f </em>= <em>ρ</em><em>u</em><strong>¨                                                              </strong>(1)

can be put in its weak form. By redefining the second derivative and converting to index notation, the governing equation becomes

(2)

which can be rearranged to

(3)

By integrating over the body Ω, the integral becomes

(4)

We can multiply by a perturbation <em>v<sub>i </sub></em>and apply the chain rule to the first term.

(5)

Applying the divergence theorem yields

(6)

The governing equation in its weak form becomes

Length terms, traction and stress terms, and time terms in the weak form equation can be nondimensionalized by the beam length <em>L</em>, Young’s modulus <em>E</em>, and characteristic time <em>t</em>˜.

<h1>Part 1b: Just Displacement</h1>

Deflection was determined at the center of a beam with length <em>L </em>= 1m, width <em>W </em>= 0<em>.</em>2m, height <em>H </em>= 0<em>.</em>2m, Youngs Modulus <em>E </em>= 200GPa, Poisson ratio <em>ν </em>= 0<em>.</em>3, and density 7800kg m<sup>−3</sup>. The beam is clamped at both ends, and a downward force of <em>F </em>= 100N is applied over the area on the top surface near the center, within 0.2 m along the width, and within 0.02 m along the length of the beam.

The deflection was determined to be 3<em>.</em>75 × 10<sup>−8 </sup>m downward, and the Paraview output is shown in Fig 1

Figure 1: Paraview output for problem 1b

The code used for part 1b can be seen in <strong>Appendix 1b</strong>.

<h1>Part 1c: Free Vibration</h1>

The force was removed, and the calculated deflection was set as an initial condition of the beam. A free vibration was modeled, and a time-plot can be seen in Fig 2.

The period of oscillation was determined to be 0.0011 seconds, and the natural frequency was determined to be 5681 rad/s or 904.2 Hz.

Figure 2: Vertical Deflection vs Time

The code used for Problem 1c is shown in <strong>Appendix 1c</strong>.

<h1>Part 1d: Natural Frequency, Changing Dimensions</h1>

<strong>Changing Width</strong>

Next, the beam width was varied to assess how natural frequency changes with beam width. Widths of 0.2, 0.4, 0.6, 0.8, and 1.0 m were assessed. A plot of natural frequencies over beam width is shown in Fig 3.

The Python code used for part 1d-i is shown in <strong>Appendix 1d-i</strong>.

Figure 3: Natural Frequency over Beam Width for Problem 1d

<strong>Changing Height</strong>

Next, the beam height was varied to assess how natural frequency changes with beam height. Heights of 0.2, 0.4, 0.6, 0.8, and 1.0 m were assessed. A plot of natural frequencies over beam Height is shown in Fig 4.

The Python code used for part 1d-ii is shown in <strong>Appendix 1d-ii</strong>.

Figure 4: Natural Frequency over Beam Height for Problem 1d

<h1>Appendix 1b: Code for Problem 1b</h1>

<em>”””</em>

<em>Python        script      for      Part 1b     of     Project     2a</em>

<em>Original          Author :         Vinamra Agrawal</em>

<em>Date :                            January    25 ,    2019</em>

<em>Edited By:                            Omkar Mulekar</em>

<em>Date :                            February     10 ,    2019</em>

<em>”””</em>

<strong>from   </strong>future          <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>import </strong>matplotlib

matplotlib . use (”Agg”)

<strong>import </strong>matplotlib . pyplot as plt <strong>from </strong>ufl <strong>import </strong>nabla div <strong>import </strong>math

<em>#============================================================== # Define       System           Properties</em>

<em>#============================================================== </em>length = 1;

W = 0.2;

H = 0.2;

a = 0.04/ length ; b = 0.4∗H/length ;

area = a∗b; F = −100

youngs = 200e9 <em># Youngs </em>nu = 0.3 <em># Poisson </em>rho = 7800 <em># Density</em>

<em># Lame parameters </em>mu = (youngs)/(2∗(1+nu) )

lambda          = (nu∗youngs)/((1+nu)∗(1−2∗nu) )

g = 10

traction applied = F/area

<em>#============================================================== #      Dimensionless           parameters</em>

<em>#============================================================== </em>l nd = length/length

w nd = W/length h nd = H/length

mu nd = mu/youngs lambda nd = lambda /youngs traction nd = traction applied /youngs

<em>#============================================================ # Boundaries and Geometry</em>

<em>#============================================================ </em>mesh = BoxMesh( Point (0 ,0 ,0) , Point ( l nd , w nd , h nd) ,20 ,6 ,6) V = VectorFunctionSpace (mesh , ’P’ ,1) tol = 1E−14

<strong>def </strong>boundary left (x , on boundary) : <strong>return </strong>(onboundary <strong>and </strong>near (x [0] ,0 , tol ) )

<strong>def </strong>boundary right (x , on boundary) : <strong>return </strong>on boundary <strong>and </strong>near (x [0] , l nd , tol )

bcleft = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundaryleft ) bc right = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundary right )

<em>#============================================================ </em><strong>def </strong>epsilon (u) : <strong>return </strong>0.5∗( nabla grad (u) + nabla grad (u) .T)

<strong>def </strong>sigma(u) :

<strong>return </strong>lambda nd∗nabla div (u)∗Identity (d) + mu nd∗( epsilon (u) + epsilon (u) .T)

<em>#============================================================ # First we solve the problem of a cantelever beam under fixed # load .</em>

<em>#============================================================</em>

u init = TrialFunction (V)

d = u init . geometric dimension () v = TestFunction (V)

f = Constant ((0.0 ,0.0 ,0.0) )

T init = Expression (( ’ 0.0 ’ , ’x [0] <em>&gt;</em>= 0.48∗ l &amp;&amp; x [0] <em>&lt;</em>= .52∗ l &amp;&amp; near (x [1] ,w) &amp;&amp; x [2] <em>&gt;</em>= 0.3∗h &amp;&amp; x [2] <em>&lt;</em>= 0.7∗h? A : 0.0 ’ , ’ 0.0 ’

) , degree=1, l=l nd , w=w nd ,h=h nd , A=traction nd )

F init = inner (sigma( u init ) , epsilon (v) )∗dx − dot( f , v)∗dx − dot(

T init , v)∗ds a init , L init = lhs ( F init ) , rhs ( F init )

<strong>print</strong>(”Solving the    i n i t i a l         cantelever problem”) u init = Function (V)

solve ( a init==L init , u init , [ bc left , bc right ])

w nd = u init ( lnd /2.0 ,w nd/2.0 , h nd /2.0)

w = w nd ∗ length <strong>print</strong>(w[1])

vtkfile u = File ( ’ deflection . pvd ’ ) vtkfile u <em>&lt;&lt; </em>u init

<h1>Appendix 1c: Code for Problem 1c</h1>

<em>”””</em>

<em>Python        script      for      Part 1c     of     Project     2a</em>

<em>Original          Author :         Vinamra Agrawal</em>

<em>Date :                            January    25 ,    2019</em>

<em>Edited By:                            Omkar Mulekar</em>

<em>Date :                            February     28 ,    2019</em>

<em>”””</em>

<strong>from   </strong>future          <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>import </strong>matplotlib

matplotlib . use (”Agg”)

<strong>import </strong>matplotlib . pyplot as plt <strong>from </strong>ufl <strong>import </strong>nabla div <strong>import </strong>math <strong>import </strong>numpy as np <strong>from </strong>scipy . signal <strong>import </strong>argrelextrema

<em>#============================================================== # Define       System           Properties</em>

<em>#============================================================== </em>length = 1;

W = 0.2;

H = 0.2;

a = 0.04∗ length ; b = 0.4∗H; area = a∗b; F = −100

youngs = 200e9 <em># Youngs </em>nu = 0.3 <em># Poisson </em>rho = 7800 <em># Density</em>

<em># Lame parameters </em>mu = (youngs)/(2∗(1+nu) )

lambda          = (nu∗youngs)/((1+nu)∗(1−2∗nu) )

g = 10

traction applied = F/area

<em>#============================================================== #      Dimensionless           parameters</em>

<em>#============================================================== </em>lnd = length/length

wnd = W/length hnd = H/length

bar speed = math. sqrt (youngs/rho) t char = length/bar speed t = 0 t i = 0.5

dt = 0.1 num steps = 150

mu nd = mu/youngs lambda nd = lambda /youngs traction nd = traction applied /youngs

<em>#============================================================ # Boundaries and Geometry</em>

<em>#============================================================ </em>mesh = BoxMesh( Point (0 ,0 ,0) , Point ( l nd , w nd , h nd) ,20 ,6 ,6) V = VectorFunctionSpace (mesh , ’P’ ,1) tol = 1E−14

<strong>def </strong>boundary left (x , on boundary) : <strong>return </strong>(onboundary <strong>and </strong>near (x [0] ,0 , tol ) )

<strong>def </strong>boundary right (x , on boundary) : <strong>return </strong>on boundary <strong>and </strong>near (x [0] , l nd , tol )

bcleft = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundaryleft ) bc right = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundary right )

<em>#============================================================ </em><strong>def </strong>epsilon (u) :

<strong>return </strong>0.5∗( nabla grad (u) + nabla grad (u) .T)

<strong>def </strong>sigma(u) :

<strong>return </strong>lambda nd∗nabla div (u)∗Identity (d) + mu nd∗( epsilon (u) + epsilon (u) .T)

<em>#============================================================ # First we solve the problem of a cantelever beam under fixed # load .</em>

<em>#============================================================ </em>u init = TrialFunction (V)

d = u init . geometric dimension () v = TestFunction (V)

f = Constant ((0.0 ,0.0 ,0.0) )

T init = Expression (( ’ 0.0 ’ , ’x [0] <em>&gt;</em>= 0.48∗ l &amp;&amp; x [0] <em>&lt;</em>= .52∗ l &amp;&amp; near (x [1] ,w) &amp;&amp; x [2] <em>&gt;</em>= 0.3∗h &amp;&amp; x [2] <em>&lt;</em>= 0.7∗h? A : 0.0 ’ , ’ 0.0 ’

) , degree=1, l=l nd , w=w nd ,h=h nd , A=traction nd )

F init = inner (sigma( u init ) , epsilon (v) )∗dx − dot( f , v)∗dx − dot(

T init , v)∗ds a init , L init = lhs ( F init ) , rhs ( F init )

<strong>print</strong>(”Solving the i n i t i a l cantelever problem”) u init = Function (V) solve ( a init==L init , u init , [ bc left , bc right ])

w nd = u init ( lnd /2.0 ,w nd/2.0 , h nd /2.0)

w = w nd ∗ length <strong>print</strong>(w[1])

<em>#============================================================</em>

<em># Next we use            this     as   i n i t i a l      condition ,      l e t    the     force      go and</em>

<em># study        the      vertical         vibrations      of     the beam</em>

<em>#============================================================ </em>u n = interpolate (Constant ((0.0 ,0.0 ,0.0) ) ,V) u n 1 = interpolate (Constant ((0.0 ,0.0 ,0.0) ) ,V) u n . assign ( u init ) u n 1 . assign ( u init ) T n = Constant ((0.0 ,0.0 ,0.0) )

u = TrialFunction (V) d = u. geometric dimension () v = TestFunction (V)

F = (dt∗dt)∗inner (sigma(u) , epsilon (v) )∗dx + dot(u, v)∗dx − (dt∗dt)∗ dot( f , v)∗dx − (dt∗dt)∗dot (T n , v)∗ds − 2.0∗ dot(u n , v)∗dx + dot( u n 1 , v)∗dx

a ,L = lhs (F) , rhs (F)

xdmffileu = XDMFFile( ’ results / solution . xdmf ’ ) xdmffiles = XDMFFile( ’ results / stress . xdmf ’ ) u = Function (V)

u store = [0] ∗ num steps time = [0] ∗ num steps

index = 0 <strong>for </strong>n <strong>in range</strong>( num steps ) : <strong>print</strong>(”time = %.2f ” % t ) Tn . t = t

solve (a == L, u, [ bc left , bc right ]) ugrab = u(0.5 ,0.1 ,0.1) u store [n] = u grab [1]

<strong>if </strong>(<strong>abs</strong>(t−index ) <em>&lt;</em>0.01) :

<strong>print</strong>(”Writing output          f i l e s . . . ”)

xdmffile u . write (u∗length , t )

W = TensorFunctionSpace (mesh , ”Lagrange” , 1) stress = lambda ∗nabla div (u)∗Identity (d) + mu∗( epsilon (u)

+ epsilon (u) .T) xdmffile s . write ( project ( stress ,W) , t ) index += 1

time [n] = t t+=dt

u n 1 . assign (u n) u n . assign (u)

<em># Get period of o s c i l l a t i o n </em>u np = np. array ( u store ) min args = argrelextrema (u np ,np. less ) period = ( time [ min args [ 0 ] [ 1 ] ] − time [ min args [ 0 ] [ 0 ] ] ) ∗t char nat freq = 2∗math. pi /period

<strong>print</strong>(”Period of Oscillation ” , period , ” seconds”) <strong>print</strong>(”Natural Frequency : ” , nat freq , ” rad/s”)

plt . figure (1) plt . plot (time , u store ) plt . xlabel ( ’ time           [ s ] ’ ) plt . ylabel ( ’ Vertical Deflection       [m] ’ ) plt . savefig ( ’1 cfig . png ’ )

<h1>Appendix 1d-i: Code for Problem 1d-i</h1>

<em>”””</em>

<em>Python        script      for      Part 1d . i      of     Project     2a</em>

<em>Original          Author :         Vinamra Agrawal</em>

<em>Date :                            January    25 ,    2019</em>

<em>Edited By:                            Omkar Mulekar</em>

<em>Date :                            February     28 ,    2019</em>

<em>”””</em>

<strong>from   </strong>future          <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>import </strong>matplotlib

matplotlib . use (”Agg”)

<strong>import </strong>matplotlib . pyplot as plt <strong>from </strong>ufl <strong>import </strong>nabla div <strong>import </strong>math <strong>import </strong>numpy as np <strong>from </strong>scipy . signal <strong>import </strong>argrelextrema

<em>#============================================================== # Define       System           Properties</em>

<em>#============================================================== </em>length = 1;

W1 = 0.2; H = 0.2;

alpha = np. linspace (0.2 ,1 ,num=5) <strong>print</strong>(”alpha = ” , alpha )

W = alpha ∗ W1

a = 0.04∗ length ; b = 0.4∗H; area = a∗b; F = −100

youngs = 200e9 <em># Youngs </em>nu = 0.3 <em># Poisson </em>rho = 7800 <em># Density</em>

<em># Lame parameters </em>mu = (youngs)/(2∗(1+nu) )

lambda          = (nu∗youngs)/((1+nu)∗(1−2∗nu) )

g = 10

traction applied = F/area

nat freq = [0] ∗ <strong>len</strong>( alpha )

<strong>print</strong>(”Beginning Loop . . . ”) <strong>for </strong>i <strong>in range</strong>(<strong>len</strong>( alpha ) ) : <strong>print</strong>(”alpha = ” , alpha [ i ])

<em>#============================================================== #     Dimensionless parameters</em>

<em>#============================================================== </em>lnd = length/length

wnd = W[ i ]/ length hnd = H/length

bar speed = math. sqrt (youngs/rho) t char = length/bar speed t = 0 t i = 0.5

dt = 0.1 num steps = 275

mu nd = mu/youngs lambda nd = lambda /youngs traction nd = traction applied /youngs

<em>#============================================================ # Boundaries and Geometry</em>

<em>#============================================================ </em>mesh = BoxMesh( Point (0 ,0 ,0) , Point ( l nd , w nd , h nd) ,20 ,6 ,6) V = VectorFunctionSpace (mesh , ’P’ ,1) tol = 1E−14

<strong>def </strong>boundaryleft (x , on boundary) :

<strong>return </strong>(onboundary <strong>and </strong>near (x [0] ,0 , tol ) )

<strong>def </strong>boundary right (x , on boundary) : <strong>return </strong>on boundary <strong>and </strong>near (x [0] , l nd , tol )

bcleft = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundaryleft ) bc right = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundary right )

<em>#============================================================ </em><strong>def </strong>epsilon (u) : <strong>return </strong>0.5∗( nabla grad (u) + nabla grad (u) .T)

<strong>def </strong>sigma(u) :

<strong>return </strong>lambda nd∗nabla div (u)∗Identity (d) + mu nd∗( epsilon (u

) + epsilon (u) .T)

<em>#============================================================ # First we solve the problem of a cantelever beam under fixed # load .</em>

<em>#============================================================ </em>u init = TrialFunction (V)

d = u init . geometric dimension () v = TestFunction (V)

f = Constant ((0.0 ,0.0 ,0.0) )

T init = Expression (( ’ 0.0 ’ , ’x [0] <em>&gt;</em>= 0.48∗ l &amp;&amp; x [0] <em>&lt;</em>= .52∗ l &amp;&amp; near (x [1] ,w) &amp;&amp; x [2] <em>&gt;</em>= 0.3∗h &amp;&amp; x [2] <em>&lt;</em>= 0.7∗h? A : 0.0 ’ , ’

0.0 ’ ) , degree=1, l=l nd , w=w nd ,h=h nd , A=traction nd )

F init = inner (sigma( uinit ) , epsilon (v) )∗dx − dot( f , v)∗dx − dot

( T init , v)∗ds a init , L init = lhs ( F init ) , rhs ( F init )

<strong>print</strong>(”Solving the i n i t i a l cantelever problem”) u init = Function (V) solve ( a init==L init , u init , [ bc left , bc right ])

down nd = u init ( l nd /2.0 ,w nd/2.0 , h nd /2.0)

w = down nd ∗ length

<strong>print</strong>(” I n i t i a l          Displacement : ” ,w[1])

<em>#============================================================</em>

<em># Next we use        this     as   i n i t i a l      condition ,      l e t    the     force      go and</em>

<em># study      the      vertical         vibrations      of     the beam</em>

<em>#============================================================ </em>u n = interpolate (Constant ((0.0 ,0.0 ,0.0) ) ,V) u n 1 = interpolate (Constant ((0.0 ,0.0 ,0.0) ) ,V) u n . assign ( u init ) u n 1 . assign ( u init )

T n = Constant ((0.0 ,0.0 ,0.0) )

u = TrialFunction (V) d = u. geometric dimension () v = TestFunction (V)

F = (dt∗dt)∗inner (sigma(u) , epsilon (v) )∗dx + dot(u, v)∗dx − (dt∗ dt)∗dot( f , v)∗dx − (dt∗dt)∗dot (T n , v)∗ds − 2.0∗ dot(u n , v)∗dx

+ dot( u n 1 , v)∗dx a ,L = lhs (F) , rhs (F)

u store = [0] ∗ num steps time = [0] ∗ num steps

index = 0 <strong>for </strong>n <strong>in range</strong>( num steps ) : <strong>print</strong>(”time = %.2f ” % t ) T n . t = t

u = Function (V)

solve (a == L, u, [ bcleft , bc right ]) ugrab = u( l nd /2.0 ,wnd/2.0 , h nd /2.0) ustore [n] = u grab [1]

<strong>if </strong>(<strong>abs</strong>(t−index ) <em>&lt;</em>0.01) :

<em># print (” Writing        output     f i l e s . . . ” )</em>

W = TensorFunctionSpace (mesh , ”Lagrange” , 1) stress = lambda ∗nabla div (u)∗Identity (d) + mu∗( epsilon (

<ol>

 <li>u) + epsilon (u) .T)</li>

</ol>

index += 1

time [n] = t t+=dt

u n 1 . assign (u n) u n . assign (u)

plt . figure (1) plt . plot (time , u store ) plt . xlabel ( ’ time [ s ] ’ ) plt . ylabel ( ’ Vertical         Deflection       [m] ’ ) plt . savefig ( ’1 dfig test . png ’ )

<em># Get period of o s c i l l a t i o n </em>u np = np. array ( u store )




min args = argrelextrema (u np ,np. greater ) <strong>print</strong>(”min args” , min args [0]) period = ( time [ min args [ 0 ] [ 1 ] ] − time [ min args [ 0 ] [ 0 ] ] ) ∗t char nat freq [ i ] = 2∗math. pi /period

<strong>print</strong>(”Period of Oscillation ” , period , ” seconds”) <strong>print</strong>(”Natural Frequency : ” , nat freq , ” rad/s”)

plt . figure (2) plt . plot (W, nat freq , ’b−x ’ ) plt . xlabel ( ’Beam Width [m] ’ ) plt . ylabel ( ’ Natural Frequency      [ rad/s ] ’ ) plt . savefig ( ’1 dfig alpha . png ’ )

<h1>Appendix 1d-ii: Code for Problem 1d-ii</h1>

<em>”””</em>

<em>Python        script      for      Part 1d . i i      of     Project     2a</em>

<em>Original          Author :         Vinamra Agrawal</em>

<em>Date :                            January    25 ,    2019</em>

<em>Edited By:                            Omkar Mulekar</em>

<em>Date :                            February     28 ,    2019</em>

<em>”””</em>

<strong>from   </strong>future          <strong> import </strong>print function <strong>from </strong>fenics <strong>import </strong>∗ <strong>import </strong>matplotlib

matplotlib . use (”Agg”)

<strong>import </strong>matplotlib . pyplot as plt <strong>from </strong>ufl <strong>import </strong>nabla div <strong>import </strong>math <strong>import </strong>numpy as np <strong>from </strong>scipy . signal <strong>import </strong>argrelextrema

<em>#============================================================== # Define       System           Properties</em>

<em>#============================================================== </em>length = 1;

W = 0.2; H1 = 0.2; beta = np. linspace (0.2 ,1 ,num=5) <strong>print</strong>(”beta = ” , beta )

H = beta ∗ H1

<strong>print</strong>(”H = ” , H) youngs = 200e9 <em># Youngs </em>nu = 0.3 <em># Poisson </em>rho = 7800 <em># Density</em>

<em># Lame parameters </em>mu = (youngs)/(2∗(1+nu) ) lambda        = (nu∗youngs)/((1+nu)∗(1−2∗nu) )

g = 10

nat freq = [0] ∗ <strong>len</strong>( beta )

<strong>print</strong>(”Beginning Loop . . . ”) <strong>for </strong>i <strong>in range</strong>(<strong>len</strong>( beta ) ) : <strong>print</strong>(”beta = ” , beta [ i ])

<em>#============================================================== #     Dimensionless parameters</em>

<em>#============================================================== </em>lnd = length/length

wnd = W/length hnd = H[ i ]/ length

bar speed = math. sqrt (youngs/rho) t char = length/bar speed t = 0 t i = 0.5

dt = 0.1 num steps = 275

mu nd = mu/youngs lambda nd = lambda /youngs

F = −100 a = 0.04∗ length ; b = 0.4∗H[ i ] ;

tractionapplied = F/(a∗b) traction nd          = traction applied /youngs

<em>#============================================================ # Boundaries and Geometry</em>

<em>#============================================================ </em>mesh = BoxMesh( Point (0 ,0 ,0) , Point ( l nd , w nd , h nd) ,20 ,6 ,6) V = VectorFunctionSpace (mesh , ’P’ ,1) tol = 1E−14

<strong>def </strong>boundaryleft (x , on boundary) :

<strong>return </strong>(onboundary <strong>and </strong>near (x [0] ,0 , tol ) )

<strong>def </strong>boundary right (x , on boundary) : <strong>return </strong>on boundary <strong>and </strong>near (x [0] , l nd , tol )

bcleft = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundaryleft ) bc right = DirichletBC (V, Constant ((0 ,0 ,0) ) , boundary right )

<em>#============================================================ </em><strong>def </strong>epsilon (u) : <strong>return </strong>0.5∗( nabla grad (u) + nabla grad (u) .T)

<strong>def </strong>sigma(u) :

<strong>return </strong>lambda nd∗nabla div (u)∗Identity (d) + mu nd∗( epsilon (u

) + epsilon (u) .T)

<em>#============================================================ # First we solve the problem of a cantelever beam under fixed # load .</em>

<em>#============================================================ </em>u init = TrialFunction (V)

d = u init . geometric dimension () v = TestFunction (V)

f = Constant ((0.0 ,0.0 ,0.0) )

T init = Expression (( ’ 0.0 ’ , ’x [0] <em>&gt;</em>= 0.48∗ l &amp;&amp; x [0] <em>&lt;</em>= .52∗ l &amp;&amp; near (x [1] ,w) &amp;&amp; x [2] <em>&gt;</em>= 0.3∗h &amp;&amp; x [2] <em>&lt;</em>= 0.7∗h? A : 0.0 ’ , ’

0.0 ’ ) , degree=1, l=l nd , w=w nd ,h=h nd , A=traction nd )

F init = inner (sigma( uinit ) , epsilon (v) )∗dx − dot( f , v)∗dx − dot

( T init , v)∗ds a init , L init = lhs ( F init ) , rhs ( F init )

<strong>print</strong>(”Solving the i n i t i a l cantelever problem”) u init = Function (V) solve ( a init==L init , u init , [ bc left , bc right ])

down nd = u init ( l nd /2.0 ,w nd/2.0 , h nd /2.0)

w = down nd ∗ length

<strong>print</strong>(” I n i t i a l          Displacement : ” ,w[1])

<em>#============================================================</em>

<em># Next we use        this     as   i n i t i a l      condition ,      l e t    the     force      go and</em>

<em># study      the      vertical         vibrations      of     the beam</em>

<em>#============================================================ </em>u n = interpolate (Constant ((0.0 ,0.0 ,0.0) ) ,V) u n 1 = interpolate (Constant ((0.0 ,0.0 ,0.0) ) ,V) u n . assign ( u init )

un 1 . assign ( u init ) Tn = Constant ((0.0 ,0.0 ,0.0) )

u = TrialFunction (V) d = u. geometric dimension () v = TestFunction (V)

F = (dt∗dt)∗inner (sigma(u) , epsilon (v) )∗dx + dot(u, v)∗dx − (dt∗ dt)∗dot( f , v)∗dx − (dt∗dt)∗dot (T n , v)∗ds − 2.0∗ dot(u n , v)∗dx

+ dot( u n 1 , v)∗dx a ,L = lhs (F) , rhs (F)

u store = [0] ∗ num steps time = [0] ∗ num steps

index = 0 <strong>for </strong>n <strong>in range</strong>( num steps ) : <strong>print</strong>(”time = %.2f ” % t ) T n . t = t

u = Function (V)

solve (a == L, u, [ bcleft , bc right ]) ugrab = u( l nd /2.0 ,wnd/2.0 , h nd /2.0) ustore [n] = u grab [1]

<strong>if </strong>(<strong>abs</strong>(t−index ) <em>&lt;</em>0.01) :

W = TensorFunctionSpace (mesh , ”Lagrange” , 1) stress = lambda ∗nabla div (u)∗Identity (d) + mu∗( epsilon (

<ol>

 <li>u) + epsilon (u) .T)</li>

</ol>

index += 1

time [n] = t t+=dt

u n 1 . assign (u n) u n . assign (u)

plt . figure (1) plt . plot (time , u store ) plt . xlabel ( ’ time [ s ] ’ ) plt . ylabel ( ’ Vertical            Deflection       [m] ’ ) plt . savefig ( ’1 dfig test2 . png ’ )

<em># Get period of o s c i l l a t i o n </em>u np = np. array ( u store )

min args = argrelextrema (u np ,np. greater ) <strong>print</strong>(”min args” , min args [0]) period = ( time [ min args [ 0 ] [ 1 ] ] − time [ min args [ 0 ] [ 0 ] ] ) ∗t char nat freq [ i ] = 2∗math. pi /period

<strong>print</strong>(”Period of Oscillation ” , period , ” seconds”) <strong>print</strong>(”Natural Frequency : ” , nat freq , ” rad/s”)

plt . figure (2) plt . plot (H, nat freq , ’b−x ’ ) plt . xlabel ( ’Beam Height            [m] ’ ) plt . ylabel ( ’ Natural Frequency      [ rad/s ] ’ ) plt . savefig ( ’1 dfig beta . png ’ )