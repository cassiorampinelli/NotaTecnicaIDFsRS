library(shiny)
library(ggplot2)
library(sf)
library(dplyr)

library(leaflet)
library(shinythemes)
library(writexl)  # Para salvar arquivos em Excel

# Defina a interface do usuário (UI)
ui <- fluidPage(
  theme = shinytheme("flatly"),  # Tema moderno e limpo
  
  # Título do cabeçalho
  titlePanel("Curvas IDFs Interpoladas para os Municípios do RS"),
  
  # Divida a página em duas metades
  fluidRow(
    # Metade esquerda para o mapa
    column(6,
           leafletOutput("map", height = "600px")  # Aumente a altura do mapa
    ),
    
    # Metade direita para a curva IDF, tabela de dados e sobre
    column(6,
           # Seletor de município e botão de gerar curva
           selectInput("municipality", 
                       "Selecione o Município:", 
                       choices = NULL),  # As escolhas serão preenchidas pelo servidor
           actionButton("generate", "Gerar Curva IDF"),  # Botão para gerar a curva
           hr(),  # Linha horizontal para separar as seções
           tabsetPanel(
             tabPanel("Curva IDF", plotOutput("idfPlot")),  # Aba para a curva IDF
             tabPanel("Tabela de Dados", 
                      tableOutput("idfTable"),  # Aba para exibir a tabela de dados
                      downloadButton("downloadData", "Baixar Tabela de Dados")  # Botão para baixar a tabela de dados
             ),
             tabPanel("Sobre", 
                      p("Esta ferramenta permite aos usuários visualizar curvas IDF para municípios selecionados. 
                         O mapa à esquerda permite selecionar um município clicando em um ponto. 
                         Após a seleção, use o menu suspenso para gerar a curva IDF. A aba da tabela de dados 
                         fornece os dados subjacentes usados para gerar a curva IDF.")
             )  # Aba de informações ou sobre
           )
    )
  )
)

# Defina a lógica do servidor
server <- function(input, output, session) {
  
  # Carregue e prepare seus dados espaciais
  sedes.municipais_sf <- st_read("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_36_\\IDFs_RS\\sedes.municipais_sf.shp")
  
  # Converta o objeto sf para um data frame com coordenadas
  sedes_coords <- st_coordinates(sedes.municipais_sf)
  sedes.municipais_df <- cbind(sedes.municipais_sf, sedes_coords)
  
  # Atualize as opções de municípios com base nos dados
  observe({
    updateSelectInput(session, "municipality", 
                      choices = setNames(sedes.municipais_df$GEOCODIGO_, sedes.municipais_df$NAME))  # Exibe o nome do município
  })
  
  # Renderize o mapa Leaflet com pontos de municípios
  output$map <- renderLeaflet({
    leaflet(sedes.municipais_df) %>%
      addTiles() %>%
      addCircleMarkers(~X, ~Y, 
                       layerId = ~GEOCODIGO_, 
                       radius = 5, color = "blue", fillOpacity = 0.7,
                       label = ~NAME)  # Mostra o nome do município ao passar o mouse
  })
  
  # Atualize o município selecionado com base no clique no mapa
  observeEvent(input$map_marker_click, {
    selected_municipality <- input$map_marker_click$id
    updateSelectInput(session, "municipality", selected = selected_municipality)
  })
  
  # Gere a curva IDF quando o botão for clicado
  observeEvent(input$generate, {
    req(input$municipality)  # Garanta que um município esteja selecionado
    
    # Filtre o município selecionado
    selected_point <- sedes.municipais_df %>%
      filter(GEOCODIGO_ == input$municipality)
    
    # Obtenha os parâmetros IDF
    a <- selected_point$idf_a_idw
    b <- selected_point$idf_b_idw
    c <- selected_point$idf_c_idw
    d <- selected_point$idf_d_idw
    
    # Simule os dados da curva IDF
    durations <- seq(5, 24 * 60, by = 20)  # Duração em minutos
    return_periods <- c(2, 5, 25, 100, 500)  # Períodos de retorno específicos em anos
    
    idf_data <- expand.grid(Duration = durations, ReturnPeriod = return_periods)
    idf_data$Intensity <- (a * (idf_data$ReturnPeriod)^b) / ((idf_data$Duration + c)^d)
    
    # Crie o gráfico da curva IDF com fontes maiores e legenda abaixo do gráfico
    idf_plot <- ggplot(idf_data, aes(x = Duration, y = Intensity, color = as.factor(ReturnPeriod))) +
      geom_line(size = 1.2) +  # Aumente o tamanho da linha para melhor visibilidade
      labs(title = paste("Curva IDF para o Município", selected_point$NAME),
           x = "Duração (minutos)",
           y = "Intensidade (mm/h)",
           color = "Período de Retorno (anos)") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),  # Aumente o tamanho do título do eixo
        axis.text = element_text(size = 14),   # Aumente o tamanho das etiquetas do eixo
        legend.text = element_text(size = 14),  # Aumente o tamanho do texto da legenda
        legend.title = element_text(size = 16),  # Aumente o tamanho do título da legenda
        legend.position = "bottom",  # Posiciona a legenda abaixo do gráfico
        legend.box = "horizontal"  # Alinha a legenda horizontalmente
      ) +
      scale_x_continuous(
        breaks = seq(0, 24 * 60, by = 120),  # Linhas de grade principais a cada 30 minutos
        minor_breaks = seq(0, 24 * 60, by = 20)  # Linhas de grade secundárias a cada 20 minutos
      )
    
    # Exiba o gráfico
    output$idfPlot <- renderPlot({
      idf_plot
    })
    
    # Exiba a tabela de dados com nomes das colunas atualizados
    output$idfTable <- renderTable({
      idf_data %>%
        rename(
          `Duração (min)` = Duration,
          `Intensidade (mm/h)` = Intensity,
          `Período de Retorno (anos)` = ReturnPeriod
        )
    })
    
    # Função de download de dados
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("IDF_data_", selected_point$NAME, ".xlsx", sep = "")
      },
      content = function(file) {
        # Salva o arquivo com os nomes das colunas atualizados
        idf_data_renamed <- idf_data %>%
          rename(
            `Duração (min)` = Duration,
            `Intensidade (mm/h)` = Intensity,
            `Período de Retorno (anos)` = ReturnPeriod
          )
        write_xlsx(idf_data_renamed, file)
      }
    )
  })
}

# Execute a aplicação
shinyApp(ui = ui, server = server)
