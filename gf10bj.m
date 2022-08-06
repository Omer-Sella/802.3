classdef gf10bj
    % Implementation of GF(2^10) arithmetic as defined by 802.3bj generator
    % polynomial 1 + x^3 +x^10

    properties
        x0 {boolean} = 0
        x1 {boolean} = 0 
        x2 {boolean} = 0
        x3 {boolean} = 0
        x4 {boolean} = 0
        x5 {boolean} = 0
        x6 {boolean} = 0
        x7 {boolean} = 0
        x8 {boolean} = 0
        x9 {boolean} = 0
    end

    methods
        function obj = untitled(inputArg0,inputArg1,inputArg2,inputArg3,inputArg4,inputArg5,inputArg6,inputArg7,inputArg8,inputArg9)
            %Construct an instance of this class
            %Three possible options: 1. we have 10 bits 2. we have only the
            %first bit 3. No bits are given
            if nargin == 10 && islogical(inputArg0) ...
                    && islogical(inputArg1) && islogical(inputArg2)...
                    && islogical(inputArg3) && islogical(inputArg4) && ...
                    islogical(inputArg5) && islogical(inputArg6) && ...
                    islogical(inputArg7) && islogical(inputArg8) && ...
                    islogical(inputArg9)
                obj.x0 = inputArg0;
                obj.x1 = inputArg1;
                obj.x2 = inputArg2;
                obj.x3 = inputArg3;
                obj.x4 = inputArg4;
                obj.x5 = inputArg5;
                obj.x6 = inputArg6;
                obj.x7 = inputArg7;
                obj.x8 = inputArg8;
                obj.x9 = inputArg9;
            end
            if nargin == 1 && islogical(inputArg0)
                obj.x0 = inputArg0;
            end
            if (nargin ~= 0) && (nargin ~= 1) && (nargin ~= 10):
                error("Omer Sella: a gf10bj object can be initialized with either 0 (none), 1, or 10 logical values")
            end

        end
            
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end