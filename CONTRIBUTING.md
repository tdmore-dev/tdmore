# Contributing to the tdmore R package
You want to contribute to ** tdmore **  ? Great ! Thank you for your enthousiasm!

Unfortunately, we are not accepting any contributions right now. This is an academic project, supporting a PhD thesis. As soon as the thesis text is finished, this package will become open-source and (hopefully) live long and prosper.

# Technical documentation
## API for models
The API for integrating different types of models is quite straightforward. Have a look at RxODE.R and deSolve.R to find how to integrate different simulation engines.

To also integrate metadata about a model, please see nlmixr.R for an example.

## EBE routines
This will be documented in a paper.

## GGPlot integration
To know how to integrate `tdmore` with `ggplot2`, it makes sense to first know the flow for `ggplot2`. In the below, we talk about Classes and Functions: ggplot2 has a custom class system that is based on prototyping (like javascript).

- First, you build a `ggplot` object. This method is an S3 generic based on the `data` argument, and usually `ggplot.default` gets called. This method uses `fortify()` to ensure the provided data argument is a data.frame. See `fortify.lm` for an example where an `lm` object is converted by using the observed data, appended with the model predictions.
- Then, you create `geom_` and `stat_` objects. Both of these functions actually result in the same thing: a `layer()` object being created, with the layer referring to specific Stat and Geom classes.
- Then, you add the new `layer` to your `ggplot` object. This function can also be abused to add other things, e.g. a new replacement data.frame. See the `ggplot_add` S3 method.

The above structure can then be converted into graphics, either through `print(x)`, or through specialized libraries like `ggplotly()`. It is important that no special Geom, Stat or Layer objects are used; they may not be translateable by `ggplotly`.

When building the plot, the following happens:
- The data for each layer is calculated. In case the data argument is empty, the data from the ggplot is used. In case the data is a function, the data from the ggplot is filtered through that function. In case it is something else, then that data is used.
- The `Layer.setup_data` method is called. This gets access to the original plot data and all layer parameters.
- Afterwards, the plot is fully built. Scales and facets are calculated, the data.frame is severely changed, compute_aesthetics is used to fill the data.frame with `x`, `y`, ... and finally the geom's are drawn.


There are a lot of different ways to integrate ggplot2 support in this flow. Requirements:

1. Use a stable API, ideally one that is public
2. Be compatible with ggplotly and other alternate drawing mechanisms
3. Be compatible with crosstalk
4. Ideally, allow the construction `plotTemplate %+% newTdmoreFit`.

#4 is not possible with ggplot 3.0.0. With a regular StatIdentity object, we probably cannot maintain anything in the `data` argument. StatIdentity has the following interesting functions: 
- `setup_params(data, params)`
- `setup_data(data, params)`
- `compute_layer(data, scales, params)`
- `computer_group(self, data, scales)`
- `compute_panel(self, data, scales, ...)`
None of these functions get the `data` unmodified. So we cannot use `data` to pass arguments. We can only get the `plot` properties during `ggplot_add`, not during the plot_build phase.


One working method is to simply use a fake `stat_predict()` function that converts the given arguments to a call into geom_line or geom_ribbon. It can only convert into a single geom though! Furthermore, we would have to pass all tdmore arguments to every `stat_predict()` if there are multiple ones.

We can improve on this by creating a special `ggplot` object: `ggtdmore`. We do this by creating a special `ggplot.tdmorefit` function that creates a ggplot with the observations data.frame (similar to `fortify.lm`), but also adds the tdmorefit and newdata/regimen/covariates parameters to the ggtdmore object.

We can further improve. Instead of allowing stat_predict to only create a single geom, we can create a special `StatPredict` object. We can then implement `ggplot_add.StatPredict` (called during `ggplot() + stat_predict()`, where we can add multiple specific geom's to the ggplot.

Finally, we can go a little further and add a special `layer(layer_class=PredictLayer)` instead. This will allow us to override `setup_data`, so the simulation is only performed when we render the plot. This should also allow replacement of the `data` argument.
This works perfectly fine (see ggplot_integration_layer.R), but is only available in the Github version. The ggplot2 release schedule was as follows: 

| Version | Date |
------------------
| 2.0.0   | December 2015 |
| 2.1.0   | February 2016 |
| 2.2.0   | November 2016 |
| 2.2.1   | December 2016 |
| 3.0.0   | June 2018     |
| 3.1.0   | October 2018  |
| 3.2.0   | ???? |

We fear that users may wait a very long time for the release of ggplot2 3.2.0. We can depend on ggplot2 3.0.0 instead, which only offers the `ggplot_add` extension method. This means we can only perform the prediction while the plot is being defined (using `ggplot_add`), and *not* when the plot is being built. Things like `template %+% newfit` are not possible...

In conclusion:
- `ggplot.tdmorefit` converts the tdmorefit object to a data.frame using fortify: it is simply the observed data. This is passed to `ggplot`. The resulting ggplot object is modified to include all tdmorefit arguments.
- `stat_predict()` creates a special layer, with a modified setup_data that performs a new simulation.

