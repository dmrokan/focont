<div class="document">

<div class="documentwrapper">

<div class="bodywrapper">

<div class="body" data-role="main">

<div id="focont-package" class="section">

# focont package[¶](#focont-package "Permalink to this heading")

<div id="submodules" class="section">

## Submodules[¶](#submodules "Permalink to this heading")

</div>

<div id="module-focont.foc" class="section">

<span id="focont-foc-module"></span>

## focont.foc module[¶](#module-focont.foc "Permalink to this heading")

  - 
    <span class="sig-prename descclassname"><span class="pre">focont.foc.</span></span><span class="sig-name descname"><span class="pre">get\_closed\_loop\_system</span></span><span class="sig-paren">(</span>*<span class="n"><span class="pre">pdata</span></span>*,
    *<span class="n"><span class="pre">i</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">-</span>
    <span class="pre">1</span></span>*,
    *<span class="n"><span class="pre">j</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">-</span>
    <span class="pre">1</span></span>*<span class="sig-paren">)</span>[<span class="viewcode-link"><span class="pre">\[source\]</span></span>](_modules/focont/foc.html#get_closed_loop_system)[¶](#focont.foc.get_closed_loop_system "Permalink to this definition")  
    Returns the closed loop system in SciPy discrete LTI system form.
    
      - Parameters<span class="colon">:</span>
        
          - **dict** (*pdata*) – Problem data structure.
        
          - **int** (*i*) – Controller output index for the MIMO
            controller.
        
          - **int** (*j*) – Controller input index for the MIMO controller.
    
      - Return scipy.signal.lti<span class="colon">:</span>  
        SciPy (discrete) LTI system representation.
    
    Returns an m by r Python array when i or j is not provided. The ith
    row and jth column of the return value gives the discrete LTI system
    from jth input to the ith output.

<!-- end list -->

  - 
    <span class="sig-prename descclassname"><span class="pre">focont.foc.</span></span><span class="sig-name descname"><span class="pre">get\_controller</span></span><span class="sig-paren">(</span>*<span class="n"><span class="pre">pdata</span></span>*,
    *<span class="n"><span class="pre">i</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">-</span>
    <span class="pre">1</span></span>*,
    *<span class="n"><span class="pre">j</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">-</span>
    <span class="pre">1</span></span>*<span class="sig-paren">)</span>[<span class="viewcode-link"><span class="pre">\[source\]</span></span>](_modules/focont/foc.html#get_controller)[¶](#focont.foc.get_controller "Permalink to this definition")  
    Returns the controller in SciPy discrete LTI system form.
    
      - Parameters<span class="colon">:</span>
        
          - **dict** (*pdata*) – Problem data structure.
        
          - **int** (*i*) – Controller output index for the MIMO
            controller.
        
          - **int** (*j*) – Controller input index for the MIMO controller.
    
    Returns an m by r Python array when i or j is not provided. The ith
    row and jth column of the return value gives the discrete LTI system
    from jth input to the ith output.

<!-- end list -->

  - 
    <span class="sig-prename descclassname"><span class="pre">focont.foc.</span></span><span class="sig-name descname"><span class="pre">norm</span></span><span class="sig-paren">(</span>*<span class="n"><span class="pre">pdata</span></span>*,
    *<span class="n"><span class="pre">cl</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">True</span></span>*<span class="sig-paren">)</span>[<span class="viewcode-link"><span class="pre">\[source\]</span></span>](_modules/focont/foc.html#norm)[¶](#focont.foc.norm "Permalink to this definition")  
    Calculates 2-norm of the impulse response of closed or open loop
    MIMO system.
    
      - Parameters<span class="colon">:</span>
        
          - **dict** (*pdata*) – Problem data structure.
        
          - **\[TODO:type\]** (*cl*) – Calculate closed loop norm if it
            is True.
    
      - Return float<span class="colon">:</span>  
        2-norm.
    
    Proposed algorithm is supposed to minimize quadratic cost function
    of the system states.

<!-- end list -->

  - 
    <span class="sig-prename descclassname"><span class="pre">focont.foc.</span></span><span class="sig-name descname"><span class="pre">print\_results</span></span><span class="sig-paren">(</span>*<span class="n"><span class="pre">pdata</span></span>*<span class="sig-paren">)</span>[<span class="viewcode-link"><span class="pre">\[source\]</span></span>](_modules/focont/foc.html#print_results)[¶](#focont.foc.print_results "Permalink to this definition")

<!-- end list -->

  - 
    <span class="sig-prename descclassname"><span class="pre">focont.foc.</span></span><span class="sig-name descname"><span class="pre">response\_improvement</span></span><span class="sig-paren">(</span>*<span class="n"><span class="pre">pdata</span></span>*<span class="sig-paren">)</span>[<span class="viewcode-link"><span class="pre">\[source\]</span></span>](_modules/focont/foc.html#response_improvement)[¶](#focont.foc.response_improvement "Permalink to this definition")  
    Compares the 2-norms of the impulse responses of the closed loop
    system obtained by the algortihm and the open loop system if the
    open loop system is also stable.
    
      - Parameters<span class="colon">:</span>  
        **dict** (*pdata*) – Problem data structure.
    
      - Return float<span class="colon">:</span>  
        Ratio of the closed and open loop 2-norms.

<!-- end list -->

  - 
    <span class="sig-prename descclassname"><span class="pre">focont.foc.</span></span><span class="sig-name descname"><span class="pre">solve</span></span><span class="sig-paren">(</span>*<span class="n"><span class="pre">pdata</span></span>*<span class="sig-paren">)</span>[<span class="viewcode-link"><span class="pre">\[source\]</span></span>](_modules/focont/foc.html#solve)[¶](#focont.foc.solve "Permalink to this definition")  
    Solves the SOF (static output feedback) or FOC (fixed order
    controller) problem for the given LTI (discrete or continous) system
    by applying the proposed solution method \[1-2\].
    
    \[1\]: Demir, O. and Özbay, H., 2020. Static output feedback
    stabilization of discrete time linear time invariant systems based
    on approximate dynamic programming. Transactions of the Institute of
    Measurement and Control, 42(16), pp.3168-3182.
    
    \[2\]: Demir, O., 2020. Optimality based structured control of
    distributed parameter systems (Doctoral dissertation, Bilkent
    University).
    
      - Parameters<span class="colon">:</span>  
        **dict** (*pdata*) – Python dictionary of problem parameters
        obtained from
    
    system.load function of focont library.
    
    Controller is calculated by performing the following steps;
    
    1.  Find an appropriate realization of the LTI system.
    
    2\. Apply the approximate dyanmic programming (ADP) iterations to
    calculate the stabilizing controller which minimize a quadratic cost
    function similart to the well-known linear quadratic regulator (LQR)
    problem.
    
    *NOTE*: Solution is appended to the input argument pdata.

</div>

<div id="module-focont.system" class="section">

<span id="focont-system-module"></span>

## focont.system module[¶](#module-focont.system "Permalink to this heading")

  - 
    <span class="sig-prename descclassname"><span class="pre">focont.system.</span></span><span class="sig-name descname"><span class="pre">load</span></span><span class="sig-paren">(</span>*<span class="n"><span class="pre">input\_data</span></span>*<span class="sig-paren">)</span>[<span class="viewcode-link"><span class="pre">\[source\]</span></span>](_modules/focont/system.html#load)[¶](#focont.system.load "Permalink to this definition")  
    Load Fixed Order Controller problem paramters from a Python data
    structure or from a json, or mat file.
    
      - Parameters<span class="colon">:</span>  
        **dict\_or\_str** (*input\_data*) – The source from which the
        problem parameters will be loaded.
    
    input\_data can be json or mat filepath. In this case, file will be
    read and problem parameters data structre will be created from the
    json or mat file.
    
    *NOTE*: Matrices must be Python array of array of floats with
    appropriate row and column sizes (They are not numpy arrays\!). Some
    matrices can be defined as a string for ease of use. E.g:
    
    > 
    > 
    > <div>
    > 
    > C = ‘I’ or Q = ‘1e-2I’
    > 
    > They will be translated to numpy identity matrices, np.eye(n) and
    > 1e-2 \* np.eye(n), where n is the dimension of LTI systems state
    > vector.
    > 
    > </div>
    
      - focont expects the following paramters:
        
          - A: System matrix of the LTI sytem ($A \in R_{n \times
            n}$).
        
          - B: Input matrix ($B \in R_{n \times m}$).
        
          - C: Output matrix (can be defined as a string, see the
                note above.) ($C \in R_{r \times n}$).
        
          - Q (optional): Cost function weight for LTI system states
            (can be defined as a string.)
        
        Q must have the same dimension as A and must be symmetric and
        semi-positive definite. If it is not provided, its default value
        is ‘I’. 
        
          - R (optional): Cost function weight for LTI system’s
        input. (can be a string.) R must be square and have the same
        number of columns as B. It must be symmetric and positive
        definite. Its default value is ‘I’. 
          - Q0’ (optional): It has the same properties as \`Q, but it is used for calculating an appropriate realization of the LTI system as an intermediary
        step of the algorithm. Its default value is ‘I’. 
          - type (optional): It can be ‘D’ if the LTI system is discrete and ‘C’
        if it is continuous. Its default value is ‘D’. 
          - Ts (optional): It is the sampling period used for ZOH discretization of the LTI system. Its default value is ‘0.01’.
          - max\_iter (optional): Dynamic programming iterations limit. Its default value is ‘1e6’.
          - eps\_conv (optional): Condition for convergence.
        If change in the cost-to-go function is smaller than this value,
        iterations will be terminated. Its default value is ‘1e-12’.
          - zoh\_calc\_step (optional): Max number of iterations used in ZOH
        discretization. Its default value is ‘256’. 
          - structure (optional): It is ‘SOF’, if a static output feedback is wanted to be calculated. It is ‘FO’ if controller is dynamic.
          - If controller structre is dynamic, then the paramters below
            can be provided.
            
              - Ccont (optional): Output matrix of the proposed dynamic
                controller.
            
            Its default value is ‘$I_{m \times m}$’. 
              - Dcont (optional): Input to output gain of the controller. Its
            default value is ‘$0_{m \times r}$’.
              - Qcont (optional): Cost function weight on controller’s state vector. Its default value is ‘I’.
              - Q0cont (optional) 
              - Rcont (optional): Cost function weight on controller’s input
            vector. Its default value is ‘I’.

<!-- end list -->

  - 
    <span class="sig-prename descclassname"><span class="pre">focont.system.</span></span><span class="sig-name descname"><span class="pre">load\_from\_json\_file</span></span><span class="sig-paren">(</span>*<span class="n"><span class="pre">json\_filename</span></span>*<span class="sig-paren">)</span>[<span class="viewcode-link"><span class="pre">\[source\]</span></span>](_modules/focont/system.html#load_from_json_file)[¶](#focont.system.load_from_json_file "Permalink to this definition")

<!-- end list -->

  - 
    <span class="sig-prename descclassname"><span class="pre">focont.system.</span></span><span class="sig-name descname"><span class="pre">load\_from\_mat\_file</span></span><span class="sig-paren">(</span>*<span class="n"><span class="pre">filename</span></span>*<span class="sig-paren">)</span>[<span class="viewcode-link"><span class="pre">\[source\]</span></span>](_modules/focont/system.html#load_from_mat_file)[¶](#focont.system.load_from_mat_file "Permalink to this definition")

</div>

<div id="module-focont" class="section">

<span id="module-contents"></span>

## Module contents[¶](#module-focont "Permalink to this heading")

</div>

</div>

</div>

</div>

</div>

</div>

